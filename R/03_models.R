default_config <- function() list(
  tickers=c("SPY","QQQ","IWM","GLD","TLT","USO","BTC-USD","ETH-USD"),
  start=as.Date("2018-01-01"),end_exclusive=as.Date("2026-06-30"),
  validation_start=as.Date("2021-01-01"),test_start=as.Date("2022-01-01"),
  refit_every=20L,min_train=500L,df_har=8,log_floor=1e-8,
  calendar="native",seed=20261113L,version="23.0-research")
read_prices <- function(path,ticker,config) {
  x<-read.csv(path,stringsAsFactors=FALSE); x$date<-as.Date(x$date)
  cols<-c("date","open","high","low","close","adjusted","volume")
  assert(all(cols %in% names(x)),paste("Missing OHLCV columns",path))
  assert(!anyDuplicated(x$date) && !anyNA(x$date),"Duplicate/invalid price date")
  x<-x[order(x$date),];x<-x[x$date>=config$start & x$date<config$end_exclusive,]
  assert(nrow(x)>config$min_train,"Too few price observations")
  for(j in c("open","high","low","close","adjusted"))
    assert(all(is.finite(x[[j]]) & x[[j]]>0),paste(ticker,"invalid price",j))
  assert(all(x$high+1e-6>=pmax(x$open,x$close) & x$low-1e-6<=pmin(x$open,x$close)),
    paste(ticker,"OHLC ordering violation; inspect source rather than deleting"))
  x$ret<-c(NA,100*diff(log(x$adjusted)))
  x$previous_date<-c(as.Date(NA),head(x$date,-1))
  # Raw close-to-close target is NEVER clipped. Range is an intraday predictor.
  x$rs<-10000*(log(x$high/x$open)*log(x$high/x$close)+log(x$low/x$open)*log(x$low/x$close))
  assert(all(x$rs>=-1e-8),"Negative Rogers-Satchell proxy")
  x$rv<-x$ret^2; x$ticker<-ticker; x<-x[-1,];rownames(x)<-NULL;x
}
features <- function(x,floor=1e-8) {
  roll<-function(z,k) { out<-rep(NA_real_,length(z)); if(length(z)>=k)out[k:length(z)]<-vapply(k:length(z),function(i)mean(z[(i-k+1):i]),0.);out }
  lr<-log(pmax(x$rv,floor)); ls<-log(pmax(x$rs,floor))
  data.frame(d=lr,w=roll(lr,5),m=roll(lr,22),q=roll(lr,66),
    rd=ls,rw=roll(ls,5),rm=roll(ls,22))
}
har_bank <- function(x,config) {
  f<-features(x,config$log_floor); n<-nrow(x)
  target<-which(x$date>=config$validation_start)
  definitions<-list(HAR=c("d","w","m"),Range_HAR=c("rd","rw","rm"),
    HAR66=c("d","w","m","q"),Log_AR1="d")
  ans<-list(); counter<-0L
  for(model in names(definitions)) {
    cols<-definitions[[model]];xx<-as.matrix(cbind(intercept=1,f[,cols,drop=FALSE]))
    train_target<-log(pmax(x$rv,config$log_floor)); beta<-NULL; smear<-oldcorrect<-NA_real_
    last_fit<-NA_integer_; fit_end<-NA_integer_
    for(j in seq_along(target)) {
      t<-target[j];o<-t-1L
      if(j==1 || (j-1)%%config$refit_every==0) {
        usable<-which(seq_len(n)<o & complete.cases(xx))
        assert(length(usable)>=config$min_train,paste("Insufficient HAR training",model))
        fit<-lm.fit(xx[usable,,drop=FALSE],train_target[usable+1L])
        assert(fit$rank==ncol(xx) && all(is.finite(fit$coefficients)),"Singular HAR fit")
        beta<-fit$coefficients;smear<-mean(exp(fit$residuals));oldcorrect<-exp(var(fit$residuals)/2)
        last_fit<-o;fit_end<-max(usable+1L)
      }
      pred<-drop(xx[o,,drop=FALSE]%*%beta)
      v<-exp(pred)*smear
      assert(is.finite(v) && v>0,"Invalid HAR variance")
      counter<-counter+1L
      ans[[counter]]<-data.frame(ticker=x$ticker[t],origin_date=x$date[o],target_date=x$date[t],
        model=model,ret=x$ret[t],sigma2=v,df=config$df_har,
        sigma2_lognormal=exp(pred)*oldcorrect,smearing=smear,
        fit_origin=x$date[last_fit],fit_target_end=x$date[fit_end],convergence=0L)
    }
  }
  do.call(rbind,ans)
}
garch_bank <- function(x,config,models=c(GARCH="sGARCH",GJR="gjrGARCH",EGARCH="eGARCH")) {
  assert(requireNamespace("rugarch",quietly=TRUE),"Install rugarch explicitly with INSTALL_PACKAGES.R")
  assert(requireNamespace("xts",quietly=TRUE),"xts is required")
  nstart<-which(x$date>=config$validation_start)[1]-1L
  assert(nstart>=config$min_train,"Insufficient GARCH training")
  ts<-xts::xts(x$ret,order.by=x$date); out<-list()
  for(m in names(models)) {
    message(x$ticker[1]," / ",m," / refit every ",config$refit_every)
    spec<-rugarch::ugarchspec(variance.model=list(model=models[[m]],garchOrder=c(1,1)),
      mean.model=list(armaOrder=c(0,0),include.mean=FALSE),distribution.model="std")
    dir.create("data/rolling_cache",recursive=TRUE,showWarnings=FALSE)
    rollfile<-file.path("data/rolling_cache",paste0(x$ticker[1],"_",m,".rds"))
    rollkey<-list(date=x$date,returns=x$ret,nstart=nstart,refit=config$refit_every,
                  model=models[[m]],rugarch=as.character(utils::packageVersion("rugarch")))
    cached<-if(file.exists(rollfile))readRDS(rollfile) else NULL
    if(!is.null(cached)&&identical(cached$key,rollkey))z<-cached$roll else {
      z<-rugarch::ugarchroll(spec,data=ts,n.ahead=1,n.start=nstart,refit.every=config$refit_every,
        refit.window="recursive",solver="hybrid",calculate.VaR=FALSE,keep.coef=TRUE)
      saveRDS(list(key=rollkey,roll=z),rollfile)
    }
    recovery<-list(solnp=list(tol=1e-8),nlminb=list(eval.max=20000,iter.max=10000),
                   gosolnp=list(n.restarts=5,n.sim=1000,rseed=config$seed))
    for(solver in names(recovery)) if(rugarch::convergence(z)!=0) {
      message("Resuming unresolved ",x$ticker[1]," / ",m," with ",solver)
      z<-rugarch::resume(z,solver=solver,solver.control=recovery[[solver]])
      saveRDS(list(key=rollkey,roll=z),rollfile)
    }
    assert(rugarch::convergence(z)==0,paste("Unresolved GARCH windows",x$ticker[1],m))
    d<-z@forecast$density; dates<-as.Date(rownames(d))
    expected<-x$date[(nstart+1):nrow(x)]
    assert(identical(dates,expected),"GARCH date mismatch; dates are never reconstructed")
    assert(max(abs(d$Realized-x$ret[(nstart+1):nrow(x)]))<1e-10,"GARCH realized target mismatch")
    assert(all(is.finite(d$Sigma)&d$Sigma>0&is.finite(d$Shape)&d$Shape>2),"Invalid GARCH forecast")
    o<-match(dates,x$date)-1L; fitidx<-nstart+((seq_along(o)-1L)%/%config$refit_every)*config$refit_every
    out[[m]]<-data.frame(ticker=x$ticker[1],origin_date=x$date[o],target_date=dates,
      model=m,ret=d$Realized,sigma2=d$Sigma^2,df=d$Shape,sigma2_lognormal=NA_real_,
      smearing=NA_real_,fit_origin=x$date[fitidx],fit_target_end=x$date[fitidx],convergence=0L)
  }
  do.call(rbind,out)
}
check_bank <- function(b,x,config) {
  assert(!anyDuplicated(paste(b$ticker,b$model,b$target_date)),"Duplicate forecast key")
  models<-c("GARCH","GJR","EGARCH","HAR","Range_HAR","HAR66","Log_AR1")
  assert(setequal(unique(b$model),models),"Required model missing")
  for(m in models) {
    z<-b[b$model==m,];z<-z[order(z$target_date),]
    assert(identical(z$target_date,x$date[x$date>=config$validation_start]),"Incomplete forecast bank")
    assert(all(z$origin_date<z$target_date & z$fit_target_end<=z$origin_date),"Training uses future targets")
    assert(max(abs(z$ret-x$ret[match(z$target_date,x$date)]))<1e-10,"Forecast target differs from raw return")
  }
  invisible(TRUE)
}
