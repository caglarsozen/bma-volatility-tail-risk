# Same six classes and chronological selection, parameters refit every 5 instead
# of 20 native observations. Source data/targets and evaluation dates are fixed.
evaluate_revision_pool <- function(bank,cfg,full_models=c("GARCH","GJR","HAR","Range_HAR","HAR66","Log_AR1")) {
  z<-bank_matrices(bank,full_models);n<-length(z$y);m<-ncol(z$v)
  vi<-which(z$dates<cfg$test_start);ti<-which(z$dates>=cfg$test_start)
  deltas<-c(.9,.95,.97,.99,1);dms<-lapply(deltas,function(d)dma_weights(z$logf,d))
  sc<-vapply(dms,function(w)mean(mix_logscore(z$logf,w)[vi]),0.)
  hazards<-c(0,.002,.01,.05);sfs<-lapply(hazards,function(h)filter_compact(z$logf,simplex_lattice(m,5),h))
  ss<-vapply(sfs,function(s)mean(s$logscore[vi]),0.)
  sw<-static_pool(z$logf[vi,,drop=FALSE]);ww<-list(EWA=matrix(1/m,n,m),Static=matrix(rep(sw,each=n),n,m),
    DMA_097=dms[[3]],DMA_selected=dms[[which.max(sc)]],Filter=sfs[[which.max(ss)]]$weights,No_transition=sfs[[1]]$weights)
  losses<-lapply(ww,function(w){a<-evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],w[ti,,drop=FALSE]);a$date<-z$dates[ti];a$ret<-z$y[ti];a})
  p<-bank_matrices(bank,c("GARCH","GJR"));ps<-static_pool(p$logf[vi,,drop=FALSE])
  for(nm in c("Static_pair","EWA_pair")) {
    w<-if(nm=="Static_pair")ps else c(.5,.5)
    a<-evaluate_mixture(p$y[ti],p$v[ti,,drop=FALSE],p$nu[ti,,drop=FALSE],matrix(rep(w,each=length(ti)),length(ti),2));a$date<-p$dates[ti];a$ret<-p$y[ti];losses[[nm]]<-a
  }
  perf<-do.call(rbind,lapply(names(losses),function(m){a<-losses[[m]];data.frame(method=m,n=nrow(a),NLL=mean(a$NLL),QLIKE=mean(a$QLIKE),FZ0_025=mean(a$FZ_0p025))}))
  comp<-compare_losses(losses,"Filter",c("Static","DMA_selected","Static_pair"),B=999)
  list(bank=bank,losses=losses,performance=perf,comparisons=comp,selected_delta=deltas[which.max(sc)],selected_hazard=hazards[which.max(ss)],
    selection=rbind(data.frame(family="DMA",parameter=deltas,validation_LPS=sc,selected=seq_along(sc)==which.max(sc)),
      data.frame(family="Filter",parameter=hazards,validation_LPS=ss,selected=seq_along(ss)==which.max(ss))))
}
run_refit_sensitivity <- function(tickers=default_config()$tickers,outdir="results/diagnostics/refit5") {
  cfg<-default_config();cfg$refit_every<-5L;dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  root<-getwd()
  for(ticker in tickers) {
    dest<-file.path(outdir,ticker);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    if(file.exists(file.path(dest,"complete.rds"))){message("Complete refit5 cache: ",ticker);next}
    x<-read_prices(file.path("data/raw",paste0(ticker,".csv")),ticker,cfg);hb<-har_bank(x,cfg)
    setwd(dest)
    gb<-tryCatch(garch_bank(x,cfg,models=c(GARCH="sGARCH",GJR="gjrGARCH")),finally=setwd(root))
    bank<-rbind(gb,hb);ans<-evaluate_revision_pool(bank,cfg)
    saveRDS(ans,file.path(dest,"complete.rds"));write.csv(bank,file.path(dest,"forecast_bank.csv"),row.names=FALSE)
    write.csv(ans$performance,file.path(dest,"performance.csv"),row.names=FALSE)
    write.csv(ans$selection,file.path(dest,"validation_selection.csv"),row.names=FALSE)
    message("REFIT5 complete ",ticker)
  }
  writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"));invisible(outdir)
}
collect_refit_sensitivity <- function(outdir="results/diagnostics/refit5") {
  cfg<-default_config();p<-cc<-ss<-adaptation<-list()
  for(ticker in cfg$tickers) {
    a<-readRDS(file.path(outdir,ticker,"complete.rds"));p[[ticker]]<-cbind(ticker=ticker,a$performance)
    cc[[ticker]]<-cbind(ticker=ticker,a$comparisons);ss[[ticker]]<-cbind(ticker=ticker,a$selection)
    adaptation[[ticker]] <- cbind(ticker = ticker,
      compare_losses(a$losses, "Filter", "No_transition", B = 999))
  }
  c<-do.call(rbind,cc);c$p_holm_refit5_family<-p.adjust(c$p_two,"holm",n=nrow(c))
  write.csv(do.call(rbind,p),file.path(outdir,"performance.csv"),row.names=FALSE)
  write.csv(c,file.path(outdir,"comparisons.csv"),row.names=FALSE)
  write.csv(do.call(rbind,ss),file.path(outdir,"selection.csv"),row.names=FALSE)
  h <- do.call(rbind, adaptation)
  h$p_holm_refit5_transition_family <- p.adjust(h$p_two, "holm", n = nrow(h))
  write.csv(h, file.path(outdir, "adaptation_vs_no_transition.csv"), row.names = FALSE)
  invisible(outdir)
}
