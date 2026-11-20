# Controlled analytic mechanism, not a new empirical forecasting candidate.
log_proxy_constants <- function(nu=Inf) {
  if(is.infinite(nu)) {
    mu<-digamma(.5)+log(2);s2<-trigamma(.5)
  } else {
    assert(nu>2,"Unit-variance Student t needs nu>2")
    mu<-digamma(.5)-digamma(nu/2)+log(nu-2)
    s2<-trigamma(.5)+trigamma(nu/2)
  }
  inflation<-exp(mu+s2/2)
  # Gaussian conditional log score has KL equal to half the QLIKE penalty.
  qpen<-log(inflation)+1/inflation-1
  if(is.infinite(nu)) {
    kl<-qpen/2
  } else {
    s<-sqrt((nu-2)/nu);bad<-s*sqrt(inflation)
    integrand<-function(x) exp(dt(x/s,nu,log=TRUE)-log(s))*
      (dt(x/s,nu,log=TRUE)-log(s)-dt(x/bad,nu,log=TRUE)+log(bad))
    kl<-integrate(integrand,-Inf,Inf,rel.tol=1e-9,subdivisions=2000)$value
  }
  data.frame(distribution=if(is.infinite(nu))"Gaussian" else paste0("Student_t_",nu),
    nu=nu,mean_log_z2=mu,var_log_z2=s2,lognormal_variance_ratio=inflation,
    exact_smearing_factor=exp(-mu),expected_QLIKE_penalty=qpen,
    expected_logscore_penalty=kl,expected_stationary_log_odds_delta097=-kl/(1-.97),
    geometric_mean_odds_delta097=exp(-kl/(1-.97)))
}
run_log_proxy_mechanism <- function(outdir="results/log_proxy_mechanism",N=200000L) {
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  dfs<-c(Inf,2.1,2.5,3,5,8,12,30);ans<-lapply(dfs,log_proxy_constants)
  write.csv(do.call(rbind,ans),file.path(outdir,"analytic_constants.csv"),row.names=FALSE)
  # Only finite-fourth-moment cases are simulated here; this is a Monte Carlo
  # illustration with sampling uncertainty, not a PASS/FAIL numerical test.
  mc<-list();set.seed(713113)
  for(nu in c(Inf,5,8,12)) {
    z<-if(is.infinite(nu))rnorm(N) else rt(N,nu)*sqrt((nu-2)/nu)
    u<-log(z^2);m<-mean(u);v<-var(u)
    smearing<-exp(m)*mean(exp(u-m));normal<-exp(m+v/2)
    # Delta-method SE for exp(mean(log Z²)+var(log Z²)/2).
    infl<-normal*((u-m)+.5*((u-m)^2-v));se<-sd(infl)/sqrt(N)
    mc[[length(mc)+1L]]<-data.frame(distribution=if(is.infinite(nu))"Gaussian" else paste0("Student_t_",nu),
      N=N,lognormal_variance_ratio=normal,normal_MC_SE=se,
      smearing_variance_ratio=smearing,smearing_MC_SE=sd(z^2)/sqrt(N),
      smearing_identity_error=abs(smearing-mean(z^2)))
  }
  write.csv(do.call(rbind,mc),file.path(outdir,"monte_carlo_illustration.csv"),row.names=FALSE)
  invisible(outdir)
}
