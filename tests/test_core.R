run_tests <- function(outdir) {
  dir.create(outdir,showWarnings=FALSE,recursive=TRUE); checks<-list()
  test<-function(name,expr) {
    ok<-tryCatch(isTRUE(force(expr)),error=function(e)FALSE)
    checks[[length(checks)+1L]]<<-data.frame(check=name,pass=ok)
    if(!ok)stop("FAILED: ",name)
  }
  n<-300L;set.seed(987);y<-rt(n,8);v<-cbind(rep(1,n),rep(9,n));nu<-matrix(8,n,2)
  logf<-t_logpdf(y,v,nu);w<-matrix(rep(c(.7,.3),each=n),n,2)
  test("log density retains dimensions",identical(dim(logf),c(n,2L)))
  test("variance to t scale",abs(t_scale(3,8)^2*8/6-3)<1e-12)
  rr<-mix_risk(w,v,nu,.025)
  test("root residual",max(abs(rr$residual))<1e-9)
  test("ES below VaR",all(rr$es<rr$q))
  f<-function(z) .7*exp(t_logpdf(z,1,8))+.3*exp(t_logpdf(z,9,8))
  test("mixture integrates to one",abs(integrate(f,-Inf,Inf)$value-1)<1e-8)
  test("ES independent numerical integral",abs(integrate(function(z)z*f(z),-Inf,rr$q[1])$value/.025-rr$es[1])<1e-7)
  one<-mix_risk(matrix(c(1,0),1),matrix(c(1,9),1),matrix(8,1,2),.01)
  test("single component quantile",abs(one$q-t_scale(1,8)*qt(.01,8))<1e-9)
  test("alpha ordering",all(mix_risk(w,v,nu,.01)$q<rr$q & rr$q<mix_risk(w,v,nu,.05)$q))
  test("invalid weights rejected",inherits(try(mix_risk(w*2,v,nu),silent=TRUE),"try-error"))
  test("invalid df rejected",inherits(try(t_scale(1,2),silent=TRUE),"try-error"))
  test("zero return QLIKE finite",is.finite(evaluate_mixture(0,matrix(c(1,9),1),matrix(8,1,2),matrix(c(.5,.5),1))$QLIKE))
  # Compare matrix log-sum-exp to a separate scalar implementation.
  scalar<-vapply(seq_len(n),function(i)log(.7*exp(logf[i,1])+.3*exp(logf[i,2])),0.)
  test("mixture score scalar agreement",max(abs(mix_logscore(logf,w)-scalar))<1e-12)
  states<-pair_grid();sf<-state_filter(logf,states,.01)
  test("state weights sum one",max(abs(rowSums(sf$weights)-1))<1e-12)
  test("state score equals flattened mixture",max(abs(sf$logscore-mix_logscore(logf,sf$weights)))<1e-12)
  changed<-logf;changed[201:n,]<-changed[201:n,]*10
  test("state filter future invariance",max(abs(state_filter(changed,states,.01)$weights[1:201,]-sf$weights[1:201,]))==0)
  test("DMA future invariance",max(abs(dma_weights(changed)[1:201,]-dma_weights(logf)[1:201,]))==0)
  test("LDF future invariance",max(abs(ldf2_weights(changed)[1:201,]-ldf2_weights(logf)[1:201,]))==0)
  test("rolling fit future invariance",max(abs(rolling_pool(changed)[1:201,]-rolling_pool(logf)[1:201,]))==0)
  a<-build_methods(logf,1:100);b<-build_methods(changed,1:100)
  test("selection ignores evaluation targets",identical(a$selection,b$selection))
  test("static pool ignores evaluation targets",identical(a$static_weights,b$static_weights))
  for(m in names(a$weights))test(paste("all-method causality",m),max(abs(a$weights[[m]][1:201,]-b$weights[[m]][1:201,]))<1e-12)
  # Duplicate a component and split each state's component mass equally.
  dupf<-cbind(logf[,1],logf[,1],logf[,2]);dups<-cbind(states[,1]/2,states[,1]/2,states[,2])
  dup<-state_filter(dupf,dups,.01)
  test("component duplication invariance",max(abs(sf$logscore-dup$logscore))<1e-12)
  # Duplicate a state and split its prior, preserving represented hypotheses.
  prior<-rep(1/nrow(states),nrow(states));ds<-rbind(states,states[1,]);dp<-c(prior,prior[1]/2);dp[1]<-dp[1]/2
  test("state duplication invariance",max(abs(state_filter(logf,ds,.01,dp)$logscore-sf$logscore))<1e-12)
  z0<-state_filter(logf,states,0)
  test("static grid regret bound",-sum(z0$logscore)<=min(-colSums(z0$state_score))+log(nrow(states))+1e-10)
  # Exact enumeration of 3^4 HMM state paths.
  st<-states[c(1,11,21),];ll<-logf[1:4,];h<-.05;res<-state_filter(ll,st,h)
  paths<-expand.grid(rep(list(1:3),4));parts<-numeric(nrow(paths))
  for(i in seq_len(nrow(paths))) {
    path<-as.integer(paths[i,]);lp<--log(3)
    for(t in 1:4) {
      lp<-lp+res$state_score[t,path[t]]
      if(t>1)lp<-lp+log((1-h)*as.numeric(path[t]==path[t-1])+h/3)
    }
    parts[i]<-lp
  }
  test("HMM likelihood equals path enumeration",abs(sum(res$logscore)-lse(parts))<1e-12)
  test("HMM path bound",all(-sum(res$logscore)<=-parts+1e-12))
  # Censoring probabilities, including the atom outside the tail.
  cut<-rep(-1,n);cl<-censored_logf(y,logf,v,nu,cut);i<-which(y>cut)[1]
  test("censor body uses survival mass",max(abs(cl[i,]-log(c(pt(-1/t_scale(1,8),8,lower.tail=FALSE),pt(-1/t_scale(9,8),8,lower.tail=FALSE)))))<1e-12)
  test("censor tail keeps density",max(abs(cl[y<=cut,]-logf[y<=cut,]))==0)
  cv<-coverage(rep(1,100),rep(0,100),.01)
  test("zero violations flagged for independence",is.na(cv$p_ind)&&cv$status=="unidentified_transition_row")
  # Raw target and HAR future invariance with complete synthetic prices.
  nn<-1000;set.seed(15);ret<-rnorm(nn,0,1);price<-100*exp(cumsum(ret/100))
  spread<-runif(nn,.002,.025)
  px<-data.frame(date=as.Date("2018-01-01")+seq_len(nn),open=price,high=price*exp(spread),
    low=price*exp(-spread),close=price,adjusted=price,volume=1000)
  tmp<-tempfile(fileext=".csv");write.csv(px,tmp,row.names=FALSE)
  cfg<-default_config();cfg$validation_start<-as.Date("2018-01-01")+750;cfg$min_train<-500
  dat<-read_prices(tmp,"SYNTHETIC",cfg);hb<-har_bank(dat,cfg)
  px$adjusted[901:1000]<-px$adjusted[901:1000]*2;write.csv(px,tmp,row.names=FALSE)
  dat2<-read_prices(tmp,"SYNTHETIC",cfg);hb2<-har_bank(dat2,cfg);unlink(tmp)
  early<-hb$target_date<=as.Date("2018-01-01")+900
  test("HAR future price invariance",max(abs(hb$sigma2[early]-hb2$sigma2[early]))==0)
  test("raw target retains extreme return",max(dat2$ret)>60)
  test("HAR training time bound",all(hb$fit_target_end<=hb$origin_date))
  test("lognormal vs smearing non-equivalence",max(abs(hb$sigma2-hb$sigma2_lognormal))>.1)
  # KKT for a static two-component mixture, independently checked by a fine grid.
  sw<-static_pool(logf);grid<-seq(0,1,by=.001)
  val<-vapply(grid,function(x)mean(mix_logscore(logf,matrix(rep(c(x,1-x),each=n),n,2))),0.)
  test("static optimizer beats fine grid",mean(mix_logscore(logf,matrix(rep(sw,each=n),n,2)))>=max(val)-1e-9)
  out<-do.call(rbind,checks);write.csv(out,file.path(outdir,"test_results.csv"),row.names=FALSE)
  message(nrow(out)," meaningful checks PASS");invisible(out)
}
