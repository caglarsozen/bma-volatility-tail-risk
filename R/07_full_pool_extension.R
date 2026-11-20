# Exploratory extension recorded after inspecting SPY/QQQ/IWM historical results.
# Configuration is fixed here; no test-return-dependent tuning is permitted.
simplex_lattice <- function(m,units=5L) {
  compositions<-function(left,k) {
    if(k==1L)return(matrix(left,1,1))
    do.call(rbind,lapply(0:left,function(x)cbind(x,compositions(left-x,k-1L))))
  }
  compositions(as.integer(units),as.integer(m))/units
}
filter_compact <- function(logf,states,hazard) {
  check_weights(states);n<-nrow(logf);k<-nrow(states);p<-rep(1/k,k)
  w<-matrix(NA_real_,n,ncol(logf));score<-numeric(n)
  for(t in seq_len(n)) {
    w[t,]<-drop(p%*%states); mx<-max(logf[t,])
    e<-log(drop(states%*%exp(logf[t,]-mx)))+mx
    score[t]<-lse(log(p)+e);p<-(1-hazard)*softmax(log(p)+e)+hazard/k
  }
  list(weights=w,logscore=score)
}
run_full_pool_extension <- function(outdir="results/main_sample") {
  cfg<-default_config();dest<-file.path(outdir,"exploratory_full_pool")
  dir.create(dest,recursive=TRUE,showWarnings=FALSE);perf<-comp<-selection<-list()
  for(ticker in cfg$tickers) {
    b<-readRDS(file.path("data/forecast_cache",paste0(ticker,".rds")))$bank
    z<-bank_matrices(b,c("GARCH","GJR","HAR","Range_HAR","HAR66","Log_AR1"))
    vi<-which(z$dates<cfg$test_start);ti<-which(z$dates>=cfg$test_start)
    states<-simplex_lattice(6,5);hazards<-c(0,.002,.01,.05)
    candidates<-lapply(hazards,function(h)filter_compact(z$logf,states,h))
    sc<-vapply(candidates,function(x)mean(x$logscore[vi]),0.)
    w<-candidates[[which.max(sc)]]$weights
    l<-evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],w[ti,,drop=FALSE])
    l$date<-z$dates[ti];l$ret<-z$y[ti]
    # Censored full-pool variant is a secondary prespecified grid search.
    eq<-matrix(1/6,length(z$y),6);cut<-mix_risk(eq,z$v,z$nu,.10)$q
    cl<-censored_logf(z$y,z$logf,z$v,z$nu,cut)
    tails<-lapply(hazards,function(h)filter_compact(cl,states,h))
    tsc<-vapply(tails,function(x)mean(x$logscore[vi]),0.)
    tw<-tails[[which.max(tsc)]]$weights
    tl<-evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],tw[ti,,drop=FALSE])
    tl$date<-z$dates[ti];tl$ret<-z$y[ti]
    losses<-readRDS(file.path(outdir,ticker,"daily_losses.rds"))
    losses$Pool_filter_six<-l;losses$Pool_filter_six_censored<-tl
    # Same-pool static censoring and LDF are necessary target-matched benchmarks.
    cw<-static_pool(cl[vi,,drop=FALSE]);cw<-matrix(rep(cw,each=length(z$y)),length(z$y),6)
    lw<-ldf2_weights(z$logf)
    for(nm in c("Static_six_censored","LDF2_six")) {
      ww<-if(nm=="Static_six_censored")cw else lw
      ll<-evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],ww[ti,,drop=FALSE])
      ll$date<-z$dates[ti];ll$ret<-z$y[ti];losses[[nm]]<-ll
    }
    for(nm in c("Pool_filter_six","Pool_filter_six_censored","Static_six_censored","LDF2_six")) {
      ll<-losses[[nm]];perf[[length(perf)+1L]]<-data.frame(ticker=ticker,method=nm,n=nrow(ll),
        LPS=-mean(ll$NLL),NLL=mean(ll$NLL),QLIKE=mean(ll$QLIKE),FZ0_025=mean(ll$FZ_0p025))
    }
    for(candidate in c("Pool_filter_six","Pool_filter_six_censored","DMA_six","Static_six")) {
      comparators<-setdiff(c("Static_six","DMA_six","EWA_six","Static_pair","Static_six_censored","LDF2_six"),candidate)
      cc<-compare_losses(losses,candidate,comparators,B=999)
      comp[[length(comp)+1L]]<-cbind(ticker=ticker,cc)
    }
    selection[[ticker]]<-data.frame(ticker=ticker,hazard=hazards,LPS_validation=sc,
      selected=seq_along(sc)==which.max(sc),censored_validation=tsc,censored_selected=seq_along(tsc)==which.max(tsc))
    saveRDS(list(weights=w,tail_weights=tw,losses=losses),file.path(dest,paste0(ticker,".rds")))
    message("Full-pool exploratory extension: ",ticker)
  }
  cc<-do.call(rbind,comp);cc$p_holm_extension_global<-p.adjust(cc$p_two,"holm")
  write.csv(cc,file.path(dest,"comparisons.csv"),row.names=FALSE)
  write.csv(do.call(rbind,perf),file.path(dest,"performance.csv"),row.names=FALSE)
  write.csv(do.call(rbind,selection),file.path(dest,"selection.csv"),row.names=FALSE)
  invisible(dest)
}
