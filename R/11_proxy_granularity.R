# The study archive contains the derivation for this proxy-granularity diagnostic.
proxy_granularity_constants <- function(d) {
  a<-d/2;logc<-digamma(a)-log(a)+trigamma(a)/2;c<-exp(logc)
  integrand<-function(y) {f<-dnorm(y);q<-f/3+2*dnorm(y,sd=sqrt(c))/3;f*(dnorm(y,log=TRUE)-log(pmax(q,.Machine$double.xmin)))}
  k<-(1+2*c)/3
  data.frame(d=d,log_c=logc,c=c,mixture_variance=k,
    NLL_deficit=integrate(integrand,-12,12,rel.tol=1e-10)$value,
    QLIKE_deficit=log(k)+1/k-1)
}
run_proxy_granularity <- function(outdir="results/diagnostics/proxy_granularity",reps=60L,n=1600L,validation_n=400L) {
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  ds<-c(1,2,5,10,22,78,390);constants<-do.call(rbind,lapply(ds,proxy_granularity_constants))
  write.csv(constants,file.path(outdir,"analytic_granularity.csv"),row.names=FALSE)
  ans<-list();ii<-0L
  # The SAME return paths are used at every proxy granularity. Oracle population
  # retransformation constants isolate scale from information/estimation error.
  set.seed(611131L);paths<-matrix(rnorm(reps*n),n,reps)
  ti<-(validation_n+1L):n
  for(j in seq_along(ds))for(rep in seq_len(reps)) {
    y<-paths[,rep];c<-constants$c[j]
    logf<-cbind(matrix(dnorm(y,log=TRUE),n,2),matrix(dnorm(y,sd=sqrt(c),log=TRUE),n,4))
    w<-dma_weights(logf,.97);ew<-matrix(1/6,n,6)
    lnll<- -mix_logscore(logf,ew);dnll<- -mix_logscore(logf,w)
    exact<-matrix(dnorm(y,log=TRUE),n,6)
    corrected<-mix_logscore(exact,dma_weights(exact,.97))-mix_logscore(exact,ew)
    ii<-ii+1L;ans[[ii]]<-data.frame(d=ds[j],replication=rep,n=n,validation_n=validation_n,
      apparent_DMA_gain=mean(lnll[ti]-dnll[ti]),
      distorted_mass=mean(rowSums(w[ti,3:6,drop=FALSE])),
      corrected_DMA_gain=mean(corrected[ti]),corrected_max_abs_difference=max(abs(corrected[ti])))
  }
  a<-do.call(rbind,ans);write.csv(a,file.path(outdir,"paired_replications.csv"),row.names=FALSE)
  s<-do.call(rbind,lapply(ds,function(d){b<-a[a$d==d,];se<-sd(b$apparent_DMA_gain)/sqrt(nrow(b));
    data.frame(d=d,reps=nrow(b),mean_apparent_DMA_gain=mean(b$apparent_DMA_gain),MCSE=se,
      lo=mean(b$apparent_DMA_gain)-qt(.975,nrow(b)-1)*se,
      hi=mean(b$apparent_DMA_gain)+qt(.975,nrow(b)-1)*se,
      mean_distorted_mass=mean(b$distorted_mass),max_corrected_difference=max(b$corrected_max_abs_difference))}))
  write.csv(s,file.path(outdir,"simulation_summary.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
  invisible(outdir)
}
