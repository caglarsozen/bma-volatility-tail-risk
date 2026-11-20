joint_tail_wald <- function(y,q,e,alpha,lag=floor(length(y)^(1/3))) {
  g<-cbind(as.numeric(y<=q)-alpha,(e-q+(q-y)*(y<=q)/alpha)/abs(e));n<-nrow(g);m<-colMeans(g);g<-sweep(g,2,m,"-")
  V<-crossprod(g)/n
  for(j in seq_len(lag)){C<-crossprod(g[(j+1L):n,,drop=FALSE],g[1:(n-j),,drop=FALSE])/n;V<-V+(1-j/(lag+1))*(C+t(C))}
  ev<-eigen(V,symmetric=TRUE,only.values=TRUE)$values
  stat<-if(min(ev)>0)as.numeric(n*crossprod(m,solve(V,m))) else NA_real_
  data.frame(n=n,lag=lag,moment_VaR=m[1],moment_ES=m[2],wald=stat,p_two=pchisq(stat,2,lower.tail=FALSE),min_eigen=min(ev))
}
run_calibration_and_floor <- function(outdir="results/diagnostics/calibration_floor") {
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE);cfg<-default_config();base<-"results/main_sample";cal<-floors<-list()
  for(ticker in cfg$tickers) {
    l<-readRDS(file.path(base,ticker,"daily_losses.rds"));s<-readRDS(file.path(base,"reviewer_robustness",paste0(ticker,".rds")))
    l$DMA_six_selected<-s$selected_losses;l$Static_grid_six<-s$static_grid_losses
    for(m in c("Static_pair","Static_six","DMA_six","DMA_six_selected","Static_grid_six"))for(a in c(.01,.025,.05)) {
      suf<-gsub("\\.","p",as.character(a));cal[[length(cal)+1L]]<-cbind(ticker=ticker,method=m,alpha=a,
        joint_tail_wald(l$Static_pair$ret,l[[m]][[paste0("VaR_",suf)]],l[[m]][[paste0("ES_",suf)]],a))
    }
    x<-read_prices(file.path("data/raw",paste0(ticker,".csv")),ticker,cfg)
    old<-read.csv(file.path(base,ticker,"forecast_bank.csv"));for(nm in c("origin_date","target_date","fit_origin","fit_target_end"))old[[nm]]<-as.Date(old[[nm]])
    for(fl in c(1e-10,1e-8,1e-6)) {
      cf<-cfg;cf$log_floor<-fl;b<-rbind(old[old$model%in%c("GARCH","GJR"),],har_bank(x,cf))
      z<-bank_matrices(b,c("GARCH","GJR","HAR","Range_HAR","HAR66","Log_AR1"));vi<-which(z$dates<cfg$test_start);ti<-which(z$dates>=cfg$test_start);n<-length(z$y)
      sw<-static_pool(z$logf[vi,,drop=FALSE]);ds<-c(.9,.95,.97,.99,1);ww<-lapply(ds,function(d)dma_weights(z$logf,d));vs<-vapply(ww,function(w)mean(mix_logscore(z$logf,w)[vi]),0.)
      ws<-list(EWA_six=matrix(1/6,n,6),Static_six=matrix(rep(sw,each=n),n,6),DMA_six_selected=ww[[which.max(vs)]])
      for(m in names(ws)) {ll<-evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],ws[[m]][ti,,drop=FALSE],alphas=.025)
        floors[[length(floors)+1L]]<-data.frame(ticker=ticker,log_floor=fl,method=m,NLL=mean(ll$NLL),QLIKE=mean(ll$QLIKE),FZ0_025=mean(ll$FZ_0p025),selected_delta=if(m=="DMA_six_selected")ds[which.max(vs)] else NA_real_)}
    }
    message("Floor and joint calibration: ",ticker)
  }
  a<-do.call(rbind,cal);a$p_holm_joint_family<-p.adjust(a$p_two,"holm")
  write.csv(a,file.path(outdir,"joint_identification_moments.csv"),row.names=FALSE)
  write.csv(do.call(rbind,floors),file.path(outdir,"log_floor_sensitivity.csv"),row.names=FALSE)
  invisible(outdir)
}
