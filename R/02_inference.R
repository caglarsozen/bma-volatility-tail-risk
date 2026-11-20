nw_difference <- function(d,lag=NULL) {
  assert(all(is.finite(d)),"Pairwise missing losses are not silently dropped")
  n<-length(d); if(is.null(lag)) lag<-floor(n^(1/3))
  lag<-min(lag,n-1); x<-d-mean(d); lrv<-sum(x*x)/n
  if(lag>0) for(j in seq_len(lag)) lrv<-lrv+2*(1-j/(lag+1))*sum(x[(j+1):n]*x[1:(n-j)])/n
  se<-sqrt(max(0,lrv)/n)
  z<-if(se>1e-14)mean(d)/se else NA_real_
  data.frame(n=n,difference=mean(d),se=se,z=z,p_two=if(is.finite(z))2*pnorm(-abs(z)) else NA_real_,lag=lag)
}
block_interval <- function(d,B=999,block=20,seed=4815) {
  set.seed(seed); n<-length(d); k<-ceiling(n/block)
  draws<-replicate(B,{
    starts<-sample.int(n,k,replace=TRUE)
    idx<-as.vector(outer(0:(block-1),starts-1,"+"))%%n+1L
    mean(d[idx[seq_len(n)]])
  })
  unname(quantile(draws,c(.025,.975),type=8))
}
compare_losses <- function(losses,reference,competitors,metrics=c("NLL","FZ_0p025"),B=999) {
  ans<-list(); k<-0
  for(cmp in competitors) for(metric in metrics) {
    d<-losses[[reference]][[metric]]-losses[[cmp]][[metric]]
    k<-k+1; r<-nw_difference(d); ci<-block_interval(d,B=B,seed=4815+k)
    ans[[k]]<-cbind(candidate=reference,comparator=cmp,metric=metric,r,lo=ci[1],hi=ci[2])
  }
  out<-do.call(rbind,ans); out$p_holm<-p.adjust(out$p_two,"holm");out
}
coverage <- function(y,q,alpha) {
  hit<-as.integer(y<=q); n<-length(hit); x<-sum(hit)
  llb<-function(a,b,p) { if(a+b==0)return(0); dbinom(b,a+b,p,log=TRUE)-lchoose(a+b,b) }
  ll0<-llb(n-x,x,alpha); ll1<-llb(n-x,x,x/n); uc<-2*(ll1-ll0)
  a<-head(hit,-1); b<-tail(hit,-1)
  n00<-sum(a==0&b==0);n01<-sum(a==0&b==1);n10<-sum(a==1&b==0);n11<-sum(a==1&b==1)
  ind<-NA_real_; status<-"ok_asymptotic"
  if(n00+n01==0 || n10+n11==0) status<-"unidentified_transition_row" else {
    p<-(n01+n11)/(n-1)
    ind<-max(0,2*(llb(n00,n01,n01/(n00+n01))+llb(n10,n11,n11/(n10+n11))-
                       llb(n00+n10,n01+n11,p)))
  }
  ci<-binom.test(x,n,alpha)$conf.int
  data.frame(n=n,hits=x,expected=n*alpha,rate=x/n,rate_lo=ci[1],rate_hi=ci[2],
    LRuc=uc,p_uc=pchisq(uc,1,lower.tail=FALSE),LRind=ind,
    p_ind=pchisq(ind,1,lower.tail=FALSE),LRcc=uc+ind,
    p_cc=pchisq(uc+ind,2,lower.tail=FALSE),status=status,
    sparse_warning=x<20 || n11<5)
}
tail_calibration <- function(y,q,e,alpha) {
  # Joint identification moments, scaled by the magnitude of forecast ES.
  h<-as.numeric(y<=q); m1<-h-alpha
  m2<-(e-q+(q-y)*h/alpha)/abs(e)
  a<-nw_difference(m1); b<-nw_difference(m2)
  data.frame(var_moment=mean(m1),var_moment_p=a$p_two,es_moment=mean(m2),
             es_moment_p=b$p_two,interpretation="marginal HAC moment diagnostics; not an ESR test")
}
pit_diagnostics <- function(pit) {
  z<-qnorm(pmin(pmax(pit,1e-10),1-1e-10))
  ties<-sum(duplicated(pit))
  kp<-if(ties>0)NA_real_ else ks.test(pit,"punif")$p.value
  data.frame(mean_pit=mean(pit),KS_p_iid_reference=kp,PIT_ties=ties,
    KS_status=if(ties>0)"not_reported_ties" else "iid_reference_only",
    transformed_mean=mean(z),transformed_var=var(z),
    LB_z_p=Box.test(z,lag=10,type="Ljung-Box")$p.value,
    LB_z2_p=Box.test(z^2,lag=10,type="Ljung-Box")$p.value)
}
