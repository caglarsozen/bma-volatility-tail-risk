run_mechanism_audit <- function(outdir="results/main_sample") {
  cfg<-default_config();dest<-file.path(outdir,"exploratory_full_pool");summary<-comparisons<-list()
  for(ticker in cfg$tickers) {
    b<-readRDS(file.path("data/forecast_cache",paste0(ticker,".rds")))$bank
    z<-bank_matrices(b,c("GARCH","GJR","HAR","Range_HAR","HAR66","Log_AR1"))
    vi<-z$dates<cfg$test_start;ti<-!vi
    f<-filter_compact(z$logf,simplex_lattice(6,5),0)
    e<-evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],f$weights[ti,,drop=FALSE])
    obj<-readRDS(file.path(dest,paste0(ticker,".rds")));losses<-obj$losses;losses$Static_grid_six<-e
    cc<-compare_losses(losses,"Pool_filter_six",c("Static_grid_six"),B=999)
    comparisons[[ticker]]<-cbind(ticker=ticker,cc)
    summary[[ticker]]<-data.frame(ticker=ticker,method="Static_grid_six",n=nrow(e),LPS=-mean(e$NLL),
      NLL=mean(e$NLL),QLIKE=mean(e$QLIKE),FZ0_025=mean(e$FZ_0p025))
    # Prediction-weight dispersion is descriptive, not a new uncertainty index.
    ww<-static_pool(z$logf[vi,,drop=FALSE]);ff<-exp(z$logf[vi,,drop=FALSE]-apply(z$logf[vi,,drop=FALSE],1,max))
    gap<-max(colMeans(ff/drop(ff%*%ww)))-1
    write.csv(data.frame(model=colnames(z$v),static_validation_weight=ww,
      mean_filter_weight=colMeans(obj$weights[ti,,drop=FALSE]),static_duality_gap=gap),
      file.path(dest,paste0(ticker,"_weight_summary.csv")),row.names=FALSE)
  }
  cc<-do.call(rbind,comparisons)
  # Include all declared contrasts in the family, including exact forecast ties
  # whose zero-variance loss difference has no normal-test statistic.
  cc$p_holm_mechanism_family<-p.adjust(cc$p_two,"holm",n=nrow(cc))
  write.csv(cc,file.path(dest,"adaptation_vs_static_grid.csv"),row.names=FALSE)
  write.csv(do.call(rbind,summary),file.path(dest,"static_grid_performance.csv"),row.names=FALSE)
}
