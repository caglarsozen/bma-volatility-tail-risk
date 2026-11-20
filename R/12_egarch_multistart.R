# Training-likelihood-only numerical repair. No forecast cap, no target deletion.
egarch_spec <- function() rugarch::ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),
  mean.model=list(armaOrder=c(0,0),include.mean=FALSE),distribution.model="std")
egarch_fit_candidates <- function(y,previous=NULL) {
  start<-list(default=NULL,canonical=list(omega=.05*log(mean(y*y)),alpha1=-.05,beta1=.95,gamma1=.15,shape=8))
  if(!is.null(previous)) start$previous<-as.list(previous)
  fits<-list();diags<-list();k<-0L
  for(solver in c("solnp","nlminb"))for(nm in names(start)) {
    sp<-egarch_spec();if(!is.null(start[[nm]]))rugarch::setstart(sp)<-start[[nm]]
    ctl<-if(solver=="solnp")list(tol=1e-8) else list(eval.max=10000,iter.max=5000)
    ws<-character()
    fit<-withCallingHandlers(tryCatch(rugarch::ugarchfit(sp,data=y,solver=solver,solver.control=ctl),error=function(e)e),
      warning=function(w){ws<<-c(ws,conditionMessage(w));invokeRestart("muffleWarning")})
    k<-k+1L;ok<-FALSE;ll<-NA_real_;cv<-NA_integer_;co<-NULL;he<-NA_real_
    if(inherits(fit,"uGARCHfit")) {
      co<-rugarch::coef(fit);ll<-as.numeric(rugarch::likelihood(fit));cv<-rugarch::convergence(fit)
      if(length(ll)!=1L)ll<-NA_real_
      if(length(cv)!=1L)cv<-NA_integer_
      ok<-isTRUE(cv==0L)&&length(co)==5L&&all(is.finite(co))&&is.finite(ll)&&isTRUE(co["shape"]>2)&&isTRUE(abs(co["beta1"])<1)
      he<-tryCatch(min(eigen((fit@fit$hessian+t(fit@fit$hessian))/2,symmetric=TRUE,only.values=TRUE)$values),error=function(e)NA_real_)
    }
    fits[k]<-list(if(ok)fit else NULL)
    diags[[k]]<-data.frame(solver=solver,start=nm,convergence=cv,valid=ok,training_LL=ll,
      native_hessian_min=he,warning=paste(unique(ws),collapse=" | "))
  }
  d<-do.call(rbind,diags);idx<-which(d$valid)
  assert(length(idx)>0,"No admissible EGARCH training fit")
  best<-idx[which.max(d$training_LL[idx])];d$selected<-seq_len(nrow(d))==best
  list(fit=fits[[best]],diagnostics=d,coefficients=rugarch::coef(fits[[best]]))
}
run_egarch_multistart_bank <- function(tickers=default_config()$tickers,outdir="results/diagnostics/egarch_full") {
  cfg<-default_config();dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  for(ticker in tickers) {
    complete<-file.path(outdir,paste0(ticker,"_complete.rds"))
    if(file.exists(complete)){message("Completed EGARCH cache: ",ticker);next}
    x<-read_prices(file.path("data/raw",paste0(ticker,".csv")),ticker,cfg)
    first<-which(x$date>=cfg$validation_start)[1L];starts<-seq(first,nrow(x),by=cfg$refit_every)
    previous<-NULL;rows<-diags<-coeffs<-list()
    for(k in seq_along(starts)) {
      file<-file.path(outdir,paste0(ticker,"_block_",k,".rds"));o<-starts[k]-1L;last<-min(o+cfg$refit_every,nrow(x))
      if(file.exists(file))b<-readRDS(file) else {
        ans<-egarch_fit_candidates(x$ret[seq_len(o)],previous)
        fixed<-egarch_spec();rugarch::setfixed(fixed)<-as.list(ans$coefficients)
        # n.old pins initialization to the estimation sample; sigma at target t
        # uses returns through t-1. Future-mutation validation is supplied separately.
        fil<-rugarch::ugarchfilter(fixed,data=x$ret[seq_len(last)],n.old=o)
        v<-as.numeric(rugarch::sigma(fil))[(o+1L):last]^2
        assert(all(is.finite(v)&v>0),paste("Invalid repaired forecast",ticker,k))
        rr<-data.frame(ticker=ticker,origin_date=x$date[o:(last-1L)],target_date=x$date[(o+1L):last],
          model="EGARCH",ret=x$ret[(o+1L):last],sigma2=v,df=unname(ans$coefficients["shape"]),
          sigma2_lognormal=NA_real_,smearing=NA_real_,fit_origin=x$date[o],fit_target_end=x$date[o],convergence=0L)
        dd<-cbind(ticker=ticker,block=k,fit_origin=x$date[o],ans$diagnostics)
        cc<-data.frame(ticker=ticker,block=k,fit_origin=x$date[o],as.list(ans$coefficients),check.names=FALSE)
        b<-list(rows=rr,diagnostics=dd,coefficients=ans$coefficients,coefficient_row=cc,
          n_train=o,native_hessian=ans$fit@fit$hessian)
        saveRDS(b,file)
      }
      previous<-b$coefficients;rows[[k]]<-b$rows;diags[[k]]<-b$diagnostics;coeffs[[k]]<-b$coefficient_row
      if(k%%10==0||k==length(starts))message(ticker," EGARCH ",k,"/",length(starts))
    }
    ans<-list(bank=do.call(rbind,rows),diagnostics=do.call(rbind,diags),coefficients=do.call(rbind,coeffs))
    saveRDS(ans,complete);write.csv(ans$bank,file.path(outdir,paste0(ticker,"_forecast_bank.csv")),row.names=FALSE)
    write.csv(ans$diagnostics,file.path(outdir,paste0(ticker,"_training_diagnostics.csv")),row.names=FALSE)
    write.csv(ans$coefficients,file.path(outdir,paste0(ticker,"_coefficients.csv")),row.names=FALSE)
    message("EGARCH completed ",ticker," max variance ",max(ans$bank$sigma2))
  }
  writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"));invisible(outdir)
}
run_egarch_curvature_audit <- function(outdir="results/diagnostics/egarch_full",
                                      tickers=default_config()$tickers,
                                      outfile=file.path(outdir,"curvature_audit.csv")) {
  cfg<-default_config();ans<-list()
  for(ticker in tickers) {
    x<-read_prices(file.path("data/raw",paste0(ticker,".csv")),ticker,cfg)
    blocks<-readRDS(file.path(outdir,paste0(ticker,"_complete.rds")))$coefficients
    for(k in seq_len(nrow(blocks))) {
      b<-readRDS(file.path(outdir,paste0(ticker,"_block_",k,".rds")))
      he<-tryCatch(min(eigen((b$native_hessian+t(b$native_hessian))/2,symmetric=TRUE,only.values=TRUE)$values),error=function(e)NA_real_)
      vals<-c(native=he,step001=NA_real_,step0001=NA_real_)
      if(!is.finite(he)||he<=0) {
        y<-x$ret[seq_len(b$n_train)];par<-b$coefficients
        fun<-function(p) {sp<-egarch_spec();rugarch::setfixed(sp)<-as.list(setNames(p,names(par)));
          -as.numeric(rugarch::likelihood(rugarch::ugarchfilter(sp,data=y,n.old=length(y))))}
        for(j in seq_along(c(.01,.001))) {
          h<-tryCatch(numDeriv::hessian(fun,par,method.args=list(d=c(.01,.001)[j])),error=function(e)NULL)
          if(!is.null(h))vals[j+1L]<-min(eigen((h+t(h))/2,symmetric=TRUE,only.values=TRUE)$values)
        }
      }
      ans[[length(ans)+1L]]<-data.frame(ticker=ticker,block=k,fit_origin=blocks$fit_origin[k],
        native=vals[1],step001=vals[2],step0001=vals[3],positive_verified=any(vals>0,na.rm=TRUE))
    }
  }
  a<-do.call(rbind,ans);write.csv(a,outfile,row.names=FALSE)
  if(any(!a$positive_verified))warning("Unresolved EGARCH curvature windows: ",sum(!a$positive_verified))
  invisible(a)
}

# Reassess step-sensitive likelihood curvature in centered, scaled
# coordinates. A very small coefficient can otherwise receive a perturbation
# too small to resolve the likelihood's second difference.
refine_egarch_curvature <- function(outdir = "results/diagnostics/egarch_full") {
  path <- file.path(outdir, "curvature_audit.csv")
  audit <- read.csv(path)
  audit$step_sensitive <- (!is.finite(audit$native) | audit$native <= 0) &
    (!is.finite(audit$step001) | audit$step001 <= 0 |
     !is.finite(audit$step0001) | audit$step0001 <= 0)
  audit$scaled_step002 <- audit$scaled_step001 <- NA_real_
  for (j in which(audit$step_sensitive)) {
    ticker <- audit$ticker[j]
    x <- read_prices(file.path("data/raw", paste0(ticker, ".csv")), ticker, default_config())
    fit <- readRDS(file.path(outdir, paste0(ticker, "_block_", audit$block[j], ".rds")))
    y <- x$ret[seq_len(fit$n_train)]
    center <- fit$coefficients
    scales <- pmax(abs(center), .05)
    unit_objective <- function(u) {
      spec <- egarch_spec()
      rugarch::setfixed(spec) <- as.list(center + scales*u)
      -as.numeric(rugarch::likelihood(rugarch::ugarchfilter(spec, y, n.old = length(y))))
    }
    min_eigen <- function(step) {
      H <- numDeriv::hessian(unit_objective, rep(0, length(center)),
        method.args = list(eps = step, d = step))
      min(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
    }
    audit$scaled_step002[j] <- min_eigen(.002)
    audit$scaled_step001[j] <- min_eigen(.001)
  }
  audit$positive_verified <- ifelse(audit$step_sensitive,
    is.finite(audit$scaled_step002) & audit$scaled_step002 > 0 &
    is.finite(audit$scaled_step001) & audit$scaled_step001 > 0,
    audit$positive_verified)
  write.csv(audit, path, row.names = FALSE)
  invisible(audit)
}
