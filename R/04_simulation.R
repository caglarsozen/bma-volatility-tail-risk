simulate_pair <- function(scenario,n=1600,seed=1,variance_ratio=9) {
  set.seed(seed); t<-seq_len(n)
  p<-switch(scenario,stationary_mix=rep(.70,n),changing_mix=ifelse((t-1)%/%200%%2==0,.90,.20),
    stationary_single=rep(1,n),changing_winner=as.numeric((t-1)%/%200%%2==0),
    shared_misspecification=rep(.70,n),nearly_identical=rep(.70,n))
  assert(!is.null(p),"Unknown simulation scenario")
  if(scenario=="nearly_identical")variance_ratio<-1.15
  # Predictable stochastic common scale, known to both forecasting components.
  lv<-numeric(n); for(i in 2:n)lv[i]<-.97*lv[i-1]+rnorm(1,0,.08)
  v<-cbind(exp(lv),exp(lv)*variance_ratio);nu<-matrix(8,n,2)
  g<-ifelse(runif(n)<p,1L,2L)
  truth_v<-v[cbind(t,g)]; truth_df<-if(scenario=="shared_misspecification")4 else 8
  if(scenario=="shared_misspecification")truth_v<-truth_v*1.6
  y<-rt(n,truth_df)*sqrt(truth_v*(truth_df-2)/truth_df)
  list(y=y,v=v,nu=nu,true_weight=p)
}
run_simulation <- function(outdir,reps=60,n=1600,validation_n=400) {
  scenarios<-c("stationary_mix","changing_mix","stationary_single","changing_winner", "shared_misspecification","nearly_identical")
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  results<-list(); k<-0L; sel<-list()
  for(s in seq_along(scenarios)) for(r in seq_len(reps)) {
    z<-simulate_pair(scenarios[s],n,seed=100000+s*1000+r)
    logf<-t_logpdf(z$y,z$v,z$nu)
    cut<-mix_risk(matrix(.5,n,2),z$v,z$nu,.10)$q
    tail<-censored_logf(z$y,logf,z$v,z$nu,cut)
    methods<-build_methods(logf,seq_len(validation_n),tail)
    idx<-(validation_n+1):n
    losses<-lapply(methods$weights,function(w)evaluate_mixture(z$y[idx],z$v[idx,,drop=FALSE],z$nu[idx,,drop=FALSE],w[idx,,drop=FALSE],alphas=.025))
    for(m in names(losses)) {
      k<-k+1L; l<-losses[[m]]
      results[[k]]<-data.frame(scenario=scenarios[s],replicate=r,method=m,NLL=mean(l$NLL),
        FZ0=mean(l$FZ_0p025),QLIKE=mean(l$QLIKE),coverage=mean(z$y[idx]<=l$VaR_0p025),
        weight_MAE=mean(abs(methods$weights[[m]][idx,1]-z$true_weight[idx])))
    }
    if(r==1)saveRDS(list(data=z,methods=methods,losses=losses),file.path(outdir,paste0(scenarios[s],"_example.rds")))
    ss<-methods$selection;ss$scenario<-scenarios[s];ss$replicate<-r;sel[[length(sel)+1L]]<-ss
    if(r%%10==0 || r==reps)message(scenarios[s],": ",r,"/",reps)
  }
  raw<-do.call(rbind,results);write.csv(raw,file.path(outdir,"simulation_replicates.csv"),row.names=FALSE)
  write.csv(do.call(rbind,sel),file.path(outdir,"simulation_selection.csv"),row.names=FALSE)
  summary<-aggregate(cbind(NLL,FZ0,QLIKE,coverage,weight_MAE)~scenario+method,raw,mean)
  write.csv(summary,file.path(outdir,"simulation_summary.csv"),row.names=FALSE)
  comps<-list();k<-0L
  for(sc in scenarios)for(candidate in c("Pool_filter","Pool_filter_censored"))
    for(cmp in c("Static_pair","DMA_selected","LDF2","Rolling_pool","Static_censored"))for(metric in c("NLL","FZ0")) {
      a<-raw[raw$scenario==sc&raw$method==candidate,];b<-raw[raw$scenario==sc&raw$method==cmp,]
      assert(identical(a$replicate,b$replicate),"Unpaired simulation replicates")
      d<-a[[metric]]-b[[metric]];se<-sd(d)/sqrt(length(d)); crit<-qt(.975,length(d)-1)
      k<-k+1L;comps[[k]]<-data.frame(scenario=sc,candidate=candidate,comparator=cmp,metric=metric,
        reps=length(d),difference=mean(d),MC_SE=se,lo=mean(d)-crit*se,hi=mean(d)+crit*se,
        p_two=if(se>0)2*pt(-abs(mean(d)/se),df=length(d)-1) else NA_real_)
    }
  comps<-do.call(rbind,comps); comps$p_holm_global<-p.adjust(comps$p_two,"holm")
  write.csv(comps,file.path(outdir,"simulation_comparisons.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
  invisible(summary)
}
