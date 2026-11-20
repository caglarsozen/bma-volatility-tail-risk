# Additional integration checks using the actual rugarch implementation.
source("00_RUN.R")
run_garch_tests <- function() {
  set.seed(817);n<-900;vv<-numeric(n);yy<-numeric(n);vv[1]<-1
  for(t in seq_len(n)) {
    if(t>1)vv[t]<-.05+.08*yy[t-1]^2+.87*vv[t-1]
    yy[t]<-rt(1,8)*sqrt(vv[t]*6/8)
  }
  x<-data.frame(ticker="SYNTHETIC",date=as.Date("2018-01-01")+seq_len(n),ret=yy)
  cfg<-default_config();cfg$validation_start<-x$date[751];cfg$min_train<-500
  a<-garch_bank(x,cfg,models=c(GARCH="sGARCH",GJR="gjrGARCH"))
  x$ret[851:n]<-x$ret[851:n]*3
  b<-garch_bank(x,cfg,models=c(GARCH="sGARCH",GJR="gjrGARCH"))
  ix<-a$target_date<=x$date[851]
  checks<-data.frame(check=c("GARCH_future_variance_invariance","GARCH_future_df_invariance","rugarch_std_density_agreement"),
    pass=c(max(abs(a$sigma2[ix]-b$sigma2[ix]))<1e-9,
           max(abs(a$df[ix]-b$df[ix]))<1e-9,
           max(abs(t_logpdf(a$ret,a$sigma2,a$df)-log(rugarch::ddist("std",a$ret,mu=0,sigma=sqrt(a$sigma2),shape=a$df))))<1e-10))
  write.csv(checks,"results/audit/garch_checks.csv",row.names=FALSE)
  assert(all(checks$pass),"GARCH integration audit failed")
  message("3 GARCH integration checks PASS")
}
run_garch_tests()
