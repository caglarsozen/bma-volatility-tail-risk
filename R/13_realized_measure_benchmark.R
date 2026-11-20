# Separate historical SPY open-to-close experiment with observed realized kernel.
# Source scale follows Sucarrat (2021), RJ-2021-057: SPY_OC*100 and SPY_RK*100.
# Stored SPY_RK is NOT squared again. This is not the main close-to-close sample.
run_realized_measure_benchmark <- function(outdir="results/diagnostics/realized_benchmark") {
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  assert(requireNamespace("rugarch",quietly=TRUE),"rugarch is required")
  assert(requireNamespace("xts",quietly=TRUE),"xts is required to decode dataset dates")
  e<-new.env();utils::data("spyreal",package="rugarch",envir=e);d<-e$spyreal
  x<-data.frame(date=as.Date(zoo::index(d)),r=as.numeric(d[,"SPY_OC"])*100,rv=as.numeric(d[,"SPY_RK"])*100)
  assert(min(x$date)==as.Date("2002-01-02")&&max(x$date)==as.Date("2008-08-29")&&nrow(x)==1662L,"Official spyreal date/size mismatch")
  assert(all(is.finite(x$rv)&x$rv>0),"Invalid official realized measure")
  write.csv(x,file.path(outdir,"spyreal_official_scaled.csv"),row.names=FALSE)
  valstart<-as.Date("2004-01-01");teststart<-as.Date("2005-01-01");first<-which(x$date>=valstart)[1L]
  y<-xts::xts(x$r,order.by=x$date);rv<-xts::xts(x$rv,order.by=x$date)
  banks<-list()
  for(nm in c("GARCH","GJR","RealGARCH")) {
    path<-file.path(outdir,paste0(nm,"_roll.rds"))
    if(file.exists(path)) roll<-readRDS(path) else {
      mod<-c(GARCH="sGARCH",GJR="gjrGARCH",RealGARCH="realGARCH")[[nm]]
      spec<-rugarch::ugarchspec(variance.model=list(model=mod,garchOrder=c(1,1)),mean.model=list(armaOrder=c(0,0),include.mean=FALSE),distribution.model="std")
      args<-list(spec=spec,data=y,n.ahead=1,n.start=first-1L,refit.every=20,refit.window="recursive",solver="hybrid",calculate.VaR=FALSE,keep.coef=TRUE)
      if(nm=="RealGARCH")args$realizedVol<-rv
      roll<-do.call(rugarch::ugarchroll,args)
      if(rugarch::convergence(roll)!=0L)roll<-rugarch::resume(roll,solver="nlminb",solver.control=list(eval.max=20000,iter.max=10000))
      assert(rugarch::convergence(roll)==0L,paste("External fit incomplete",nm));saveRDS(roll,path)
    }
    a<-roll@forecast$density;dates<-as.Date(rownames(a));assert(identical(dates,x$date[first:nrow(x)]),"External dates")
    assert(max(abs(a$Realized-x$r[first:nrow(x)]))<1e-10,"External targets")
    banks[[nm]]<-data.frame(ticker="SPY_external",origin_date=x$date[(first-1L):(nrow(x)-1L)],target_date=dates,
      model=nm,ret=a$Realized,sigma2=a$Sigma^2,df=a$Shape,sigma2_lognormal=NA_real_)
  }
  lr<-log(x$rv);rollmean<-function(k)vapply(seq_along(lr),function(i)if(i<k)NA_real_ else mean(lr[(i-k+1L):i]),0.)
  X<-cbind(1,lr,rollmean(5),rollmean(22));har<-list()
  for(j in seq_len(nrow(x)-first+1L)) {
    t<-first+j-1L;o<-t-1L
    if(j==1L||(j-1L)%%20L==0L) {
      idx<-which(seq_len(nrow(x))<o & complete.cases(X));fit<-lm.fit(X[idx,,drop=FALSE],lr[idx+1L]);
      beta<-fit$coefficients;smear<-mean(exp(fit$residuals));old<-exp(var(fit$residuals)/2)
    }
    pred<-drop(X[o,,drop=FALSE]%*%beta)
    har[[j]]<-data.frame(ticker="SPY_external",origin_date=x$date[o],target_date=x$date[t],model="HAR_RK",ret=x$r[t],
      sigma2=exp(pred)*smear,df=8,sigma2_lognormal=exp(pred)*old)
  }
  bank<-do.call(rbind,c(banks,list(do.call(rbind,har))));write.csv(bank,file.path(outdir,"forecast_bank.csv"),row.names=FALSE)
  z<-bank_matrices(bank,c("GARCH","GJR","RealGARCH","HAR_RK"));vi<-which(z$dates<teststart);ti<-which(z$dates>=teststart);n<-nrow(z$v)
  sw<-static_pool(z$logf[vi,,drop=FALSE]);deltas<-c(.9,.95,.97,.99,1);dm<-lapply(deltas,function(d)dma_weights(z$logf,d));ds<-vapply(dm,function(w)mean(mix_logscore(z$logf,w)[vi]),0.)
  hs<-c(0,.002,.01,.05);sf<-lapply(hs,function(h)filter_compact(z$logf,simplex_lattice(4,5),h));ss<-vapply(sf,function(a)mean(a$logscore[vi]),0.)
  weights<-list(EWA=matrix(1/4,n,4),Static=matrix(rep(sw,each=n),n,4),DMA_selected=dm[[which.max(ds)]],Pool_filter=sf[[which.max(ss)]]$weights,No_transition=sf[[1]]$weights)
  for(j in 1:4){ww<-matrix(0,n,4);ww[,j]<-1;weights[[colnames(z$v)[j]]]<-ww}
  losses<-lapply(weights,function(w)evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],w[ti,,drop=FALSE]))
  oldv<-z$v;hh<-bank[bank$model=="HAR_RK",];hh<-hh[order(hh$target_date),];oldv[,4]<-hh$sigma2_lognormal
  losses$EWA_old_transform<-evaluate_mixture(z$y[ti],oldv[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],weights$EWA[ti,,drop=FALSE])
  perf<-do.call(rbind,lapply(names(losses),function(m){a<-losses[[m]];data.frame(method=m,n=nrow(a),NLL=mean(a$NLL),QLIKE=mean(a$QLIKE),FZ0_025=mean(a$FZ_0p025))}))
  cc<-do.call(rbind,lapply(c("Pool_filter","No_transition","DMA_selected"),function(m)compare_losses(losses,m,c("EWA","Static","GARCH","GJR","RealGARCH","HAR_RK"),B=999)))
  cc$p_holm_external_family<-p.adjust(cc$p_two,"holm")
  sel<-rbind(data.frame(family="DMA",parameter=deltas,validation_LPS=ds,selected=seq_along(ds)==which.max(ds)),data.frame(family="Filter",parameter=hs,validation_LPS=ss,selected=seq_along(ss)==which.max(ss)))
  write.csv(perf,file.path(outdir,"performance.csv"),row.names=FALSE);write.csv(cc,file.path(outdir,"comparisons.csv"),row.names=FALSE)
  write.csv(sel,file.path(outdir,"selection.csv"),row.names=FALSE);write.csv(data.frame(model=colnames(z$v),weight=as.numeric(sw)),file.path(outdir,"static_weights.csv"),row.names=FALSE)
  write.csv(data.frame(n=nrow(x),validation_n=length(vi),evaluation_n=length(ti),first=min(x$date),last=max(x$date),median_normal_smear_ratio=median(hh$sigma2_lognormal/hh$sigma2)),file.path(outdir,"data_audit.csv"),row.names=FALSE)
  saveRDS(list(losses=losses,weights=weights,bank=bank,dates=z$dates[ti],y=z$y[ti]),file.path(outdir,"complete.rds"))
  writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"));invisible(outdir)
}
