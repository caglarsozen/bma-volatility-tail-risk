# Predictive density pools and scoring rules.
# All means are fixed at zero; sigma2 denotes variance, not Student-t scale.
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
lse <- function(x) { m <- max(x); if (!is.finite(m)) return(m); m + log(sum(exp(x-m))) }
softmax <- function(x) exp(x-lse(x))
row_lse <- function(x) { m <- apply(x,1,max); m+log(rowSums(exp(x-m))) }
check_weights <- function(w) {
  assert(all(is.finite(w)) && all(w >= 0), "Invalid mixture weights")
  assert(max(abs(rowSums(as.matrix(w))-1)) < 1e-9, "Weights must sum to one")
}
t_scale <- function(v, nu) {
  assert(all(is.finite(v) & v > 0 & is.finite(nu) & nu > 2), "Invalid variance or df")
  sqrt(v*(nu-2)/nu)
}
t_logpdf <- function(y,v,nu) { s<-t_scale(v,nu); dt(y/s,nu,log=TRUE)-log(s) }
mix_logscore <- function(logf,w) { check_weights(w); row_lse(logf+log(w)) }
mix_cdf <- function(q,w,v,nu) {
  check_weights(w); rowSums(w*pt(q/t_scale(v,nu),nu))
}
mix_risk <- function(w,v,nu,alpha=.025,tol=1e-10) {
  w<-as.matrix(w); v<-as.matrix(v); nu<-as.matrix(nu)
  check_weights(w)
  assert(identical(dim(w),dim(v)) && identical(dim(v),dim(nu)), "Risk dimensions differ")
  assert(length(alpha)==1 && alpha>0 && alpha<.5, "Left-tail alpha required")
  sc<-t_scale(v,nu)
  # A mixture quantile is between the smallest and largest component quantiles.
  qc<-sc*qt(alpha,nu); lo<-apply(qc,1,min); hi<-apply(qc,1,max)
  for (k in 1:80) {
    mid<-(lo+hi)/2; cdf<-rowSums(w*pt(mid/sc,nu)); below<-cdf<alpha
    lo[below]<-mid[below]; hi[!below]<-mid[!below]
    if (max(hi-lo)<tol) break
  }
  q<-(lo+hi)/2; z<-q/sc
  es<-rowSums(w*(-sc*(nu+z^2)/(nu-1)*dt(z,nu)))/alpha
  resid<-rowSums(w*pt(z,nu))-alpha
  assert(max(abs(resid))<1e-8,"Mixture CDF root failed")
  assert(all(es<=q+1e-9),"ES exceeds VaR")
  list(q=q,es=es,residual=resid)
}
fz0 <- function(y,q,e,alpha=.025) {
  assert(all(e<0 & is.finite(e)), "FZ0 requires finite negative ES")
  -(y<=q)*(q-y)/(alpha*e)+q/e+log(-e)-1
}
censored_logf <- function(y,logf,v,nu,cutoff) {
  # cutoff must be known before y; same cutoff for every compared model.
  out<-logf; body<-y>cutoff
  out[body,]<-pt(cutoff[body]/t_scale(v[body,,drop=FALSE],nu[body,,drop=FALSE]),
                    nu[body,,drop=FALSE],lower.tail=FALSE,log.p=TRUE)
  out
}
static_pool <- function(logf,obs_weights=NULL,tol=1e-10,maxit=10000) {
  assert(nrow(logf)>0 && all(is.finite(logf)),"Invalid static training scores")
  m<-ncol(logf); if (m==1) return(1)
  a<-if(is.null(obs_weights)) rep(1/nrow(logf),nrow(logf)) else obs_weights/sum(obs_weights)
  f<-exp(logf-apply(logf,1,max))
  obj<-function(w) -sum(a*log(drop(f%*%w)))
  if (m==2) {
    fit<-optimize(function(x)obj(c(x,1-x)),c(0,1),tol=tol)
    x<-c(0,fit$minimum,1); val<-vapply(x,function(z)obj(c(z,1-z)),0.)
    x<-x[which.min(val)]; return(c(x,1-x))
  }
  w<-rep(1/m,m)
  for(k in seq_len(min(maxit,2000L))) {
    resp<-sweep(f,2,w,"*")/drop(f%*%w)
    wn<-colSums(resp*a)
    if(max(abs(w-wn))<tol) {w<-wn/sum(wn);break}
    w<-wn
  }
  # Correlated components can make EM weights converge extremely slowly even
  # when objective changes are tiny. Require a concave-objective duality check.
  gap<-function(x)max(colSums(f/drop(f%*%x)*a))-1
  if(gap(w)>1e-7) {
    assert(requireNamespace("nloptr",quietly=TRUE),"nloptr required for the full-pool static optimum")
    fit<-nloptr::nloptr(x0=(w+1e-8)/sum(w+1e-8),
      eval_f=function(x)list(objective=obj(x),gradient=-colSums(f/drop(f%*%x)*a)),
      lb=rep(0,m),ub=rep(1,m),eval_g_eq=function(x)sum(x)-1,
      eval_jac_g_eq=function(x)matrix(1,1,m),
      opts=list(algorithm="NLOPT_LD_SLSQP",xtol_rel=1e-11,ftol_abs=1e-13,maxeval=5000))
    candidate<-pmax(fit$solution,0);candidate<-candidate/sum(candidate)
    if(obj(candidate)<=obj(w)+1e-12)w<-candidate
  }
  assert(gap(w)<1e-6,"Static pool has unresolved optimality gap")
  attr(w,"duality_gap")<-gap(w);w/sum(w)
}
dma_weights <- function(logf,delta=.97,base=NULL,anchor=FALSE) {
  n<-nrow(logf); m<-ncol(logf); if(is.null(base)) base<-rep(1/m,m)
  assert(all(base>0) && abs(sum(base)-1)<1e-9,"DMA prior must be positive")
  # Forecast chronology: forecast w_t, then log w_{t+1}=delta log w_t+log f_t.
  # This differs from discounting the just-updated posterior; documented explicitly.
  lw<-log(base); w<-matrix(NA_real_,n,m)
  for(t in seq_len(n)) {
    w[t,]<-softmax(lw)
    lw<-delta*lw+logf[t,]+if(anchor)(1-delta)*log(base) else 0
    lw<-lw-lse(lw)
  }
  w
}
ldf2_weights <- function(logf,deltas=c(.90,.95,.97,.99,1),delta2=.99) {
  first<-lapply(deltas,function(d)dma_weights(logf,d))
  meta<-sapply(first,function(w)mix_logscore(logf,w))
  upper<-dma_weights(meta,delta2)
  out<-matrix(0,nrow(logf),ncol(logf))
  for(j in seq_along(first)) out<-out+first[[j]]*upper[,j]
  out
}
pair_grid <- function(step=.05) {
  x<-seq(0,1,by=step); assert(tail(x,1)==1,"Grid step must partition [0,1]")
  cbind(x,1-x)
}
state_filter <- function(logf,states,hazard=.01,prior=NULL) {
  # Finite-state dynamic prediction pool / fixed-share aggregation.
  # This is an adaptation of existing algorithms, not a claimed new invention.
  check_weights(states); k<-nrow(states)
  if(is.null(prior)) prior<-rep(1/k,k)
  assert(hazard>=0 && hazard<1 && all(prior>0) && abs(sum(prior)-1)<1e-9,"Invalid state prior")
  n<-nrow(logf); out<-matrix(NA_real_,n,ncol(logf)); p<-prior
  state_score<-matrix(NA_real_,n,k); logscore<-numeric(n)
  for(t in seq_len(n)) {
    out[t,]<-drop(p%*%states)
    mx<-max(logf[t,]); emission<-log(drop(states%*%exp(logf[t,]-mx)))+mx
    state_score[t,]<-emission; logscore[t]<-lse(log(p)+emission)
    post<-softmax(log(p)+emission)
    p<-(1-hazard)*post+hazard*prior
  }
  list(weights=out,logscore=logscore,state_score=state_score)
}
rolling_pool <- function(logf,window=250,refit=20) {
  n<-nrow(logf); m<-ncol(logf); w<-matrix(1/m,n,m)
  for(t in 2:n) {
    if((t-2)%%refit==0) ww<-static_pool(logf[max(1,t-window):(t-1),,drop=FALSE])
    w[t,]<-ww
  }
  w
}
build_methods <- function(logf,validation,tail_logf=NULL) {
  assert(ncol(logf)==2,"Primary identification experiment requires the fixed strong pair")
  assert(length(validation)>=30 && identical(validation,seq_len(max(validation))),"Validation must be a chronological prefix")
  n<-nrow(logf); m<-2; vc<-max(validation)
  repeat_w<-function(w) matrix(rep(w,each=n),n,m)
  sw<-static_pool(logf[validation,,drop=FALSE])
  dm_grid<-c(.90,.95,.97,.99,1)
  dms<-lapply(dm_grid,function(d)dma_weights(logf,d))
  dm_scores<-vapply(dms,function(w)mean(mix_logscore(logf,w)[validation]),0.)
  hazards<-c(0,.002,.01,.05)
  states<-pair_grid()
  sfs<-lapply(hazards,function(h)state_filter(logf,states,h))
  sf_scores<-vapply(sfs,function(z)mean(z$logscore[validation]),0.)
  methods<-list(EWA_pair=repeat_w(c(.5,.5)),Static_pair=repeat_w(sw),
    Best_single=repeat_w(as.numeric(seq_len(m)==which.max(colMeans(logf[validation,,drop=FALSE])))),
    DMA_097=dms[[3]],DMA_selected=dms[[which.max(dm_scores)]],
    Cumulative=dms[[5]],LDF2=ldf2_weights(logf),
    Rolling_pool=rolling_pool(logf),Pool_filter=sfs[[which.max(sf_scores)]]$weights,
    Static_grid=sfs[[1]]$weights)
  selection<-data.frame(family=c(rep("DMA",5),rep("Pool_filter",4)),
    parameter=c(dm_grid,hazards),validation_score=c(dm_scores,sf_scores),
    selected=c(seq_along(dm_grid)==which.max(dm_scores),seq_along(hazards)==which.max(sf_scores)))
  if(!is.null(tail_logf)) {
    ss<-static_pool(tail_logf[validation,,drop=FALSE])
    ts<-lapply(hazards,function(h)state_filter(tail_logf,states,h))
    scores<-vapply(ts,function(z)mean(z$logscore[validation]),0.)
    methods$Static_censored<-repeat_w(ss)
    methods$Pool_filter_censored<-ts[[which.max(scores)]]$weights
    selection<-rbind(selection,data.frame(family="Pool_filter_censored",parameter=hazards,
      validation_score=scores,selected=seq_along(hazards)==which.max(scores)))
  }
  # Static/selected specifications are trained with the entire validation prefix.
  # Their reported performance starts AFTER that prefix. No validation claims for them.
  list(weights=methods,selection=selection,static_weights=sw,test_start=vc+1L)
}
evaluate_mixture <- function(y,v,nu,w,alphas=c(.01,.025,.05)) {
  logf<-t_logpdf(y,v,nu); ls<-mix_logscore(logf,w)
  vv<-rowSums(w*v)
  # QLIKE in equivalent loss form log(v)+r^2/v: well defined even at r=0.
  out<-data.frame(NLL=-ls,QLIKE=log(vv)+y^2/vv,variance=vv,PIT=mix_cdf(y,w,v,nu),
    entropy=-rowSums(ifelse(w>0,w*log(pmax(w,1e-300)),0)))
  for(a in alphas) {
    z<-mix_risk(w,v,nu,a); suffix<-gsub("\\.","p",as.character(a))
    out[[paste0("VaR_",suffix)]]<-z$q; out[[paste0("ES_",suffix)]]<-z$es
    out[[paste0("FZ_",suffix)]]<-fz0(y,z$q,z$es,a)
  }
  out
}
