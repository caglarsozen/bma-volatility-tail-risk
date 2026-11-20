audit_old_means <- function(config=default_config(),outdir="results") {
  a<-lapply(config$tickers,function(ticker)read_prices(file.path("data/raw",paste0(ticker,".csv")),ticker,config))
  names(a)<-config$tickers
  common<-Reduce(intersect,lapply(a,function(x)as.character(x$date)))
  ans<-lapply(names(a),function(ticker) {
    x<-a[[ticker]];q<-quantile(x$ret,c(.001,.999),names=FALSE)
    old<-pmin(pmax(x$ret,q[1]),q[2]);old<-old[as.character(x$date)%in%common]
    data.frame(ticker=ticker,native_raw_mean=mean(x$ret),native_n=nrow(x),
      legacy_style_mean_new_vintage=mean(old),common_n=length(old),rounded4=sprintf("%.4f",mean(old)))
  })
  out<-do.call(rbind,ans);write.csv(out,file.path(outdir,"mean_return_audit.csv"),row.names=FALSE);out
}
pool_envelope <- function(y,v,nu,alphas=c(.01,.025,.05)) {
  # Every linear pool quantile lies inside the component quantile envelope.
  # Every tail-exceedance sequence is therefore bounded pathwise by these two
  # binary sequences. The envelope is a necessary, not sufficient, calibration test.
  ans<-lapply(alphas,function(alpha) {
    q<-t_scale(v,nu)*qt(alpha,nu);lo<-apply(q,1,min);hi<-apply(q,1,max)
    lower_hit<-as.integer(y<=lo);upper_hit<-as.integer(y<=hi)
    dl<-nw_difference(lower_hit-alpha);du<-nw_difference(upper_hit-alpha)
    data.frame(alpha=alpha,n=length(y),unavoidable_hits=sum(lower_hit),possible_hits=sum(upper_hit),
      unavoidable_rate=mean(lower_hit),possible_rate=mean(upper_hit),
      lower_minus_nominal=dl$difference,lower_HAC_se=dl$se,
      unavoidable_excess_p_one=if(is.finite(dl$z))pnorm(dl$z,lower.tail=FALSE) else NA_real_,
      upper_minus_nominal=du$difference,upper_HAC_se=du$se,
      insufficient_possible_hits_p_one=if(is.finite(du$z))pnorm(du$z) else NA_real_,
      exact_sample_nominal_inside=(mean(lower_hit)<=alpha&&alpha<=mean(upper_hit)))
  })
  do.call(rbind,ans)
}
run_extra_audits <- function(outdir="results/main_sample") {
  cfg<-default_config();audit_old_means(cfg)
  envelopes<-list();k<-0
  for(ticker in cfg$tickers) {
    path<-file.path(outdir,ticker,"forecast_bank.csv");if(!file.exists(path))next
    b<-read.csv(path);b$target_date<-as.Date(b$target_date);b$origin_date<-as.Date(b$origin_date)
    for(pool in c("pair","seven")) {
      z<-bank_matrices(b,if(pool=="pair")c("GARCH","GJR") else unique(b$model))
      ix<-z$dates>=cfg$test_start;k<-k+1
      envelopes[[k]]<-cbind(ticker=ticker,pool=pool,pool_envelope(z$y[ix],z$v[ix,,drop=FALSE],z$nu[ix,,drop=FALSE]))
    }
  }
  if(length(envelopes)) {
    ee<-do.call(rbind,envelopes)
    # All assets/pools/alphas/directions in one explicitly disclosed family.
    adj<-p.adjust(c(ee$unavoidable_excess_p_one,ee$insufficient_possible_hits_p_one),"holm")
    ee$unavoidable_p_holm<-head(adj,nrow(ee));ee$insufficient_p_holm<-tail(adj,nrow(ee))
    write.csv(ee,file.path(outdir,"pool_calibration_envelope.csv"),row.names=FALSE)
  }
}
