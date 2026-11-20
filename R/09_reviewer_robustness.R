# Reviewer-requested delta sensitivity and secondary QLIKE inference.
# Added after the historical experiment; no test score selects a configuration.
run_reviewer_robustness <- function(outdir="results/main_sample") {
  cfg<-default_config();dest<-file.path(outdir,"reviewer_robustness")
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  delta_tables<-qlike<-selected_comparisons<-block_tables<-calibration_tables<-repair_tables<-list()
  for(ticker in cfg$tickers) {
    bank<-readRDS(file.path("data/forecast_cache",paste0(ticker,".rds")))$bank
    z<-bank_matrices(bank,c("GARCH","GJR","HAR","Range_HAR","HAR66","Log_AR1"))
    vi<-which(z$dates<cfg$test_start);ti<-which(z$dates>=cfg$test_start)
    deltas<-c(.90,.95,.97,.99,1)
    ww<-lapply(deltas,function(d)dma_weights(z$logf,d))
    scores<-vapply(ww,function(w)mean(mix_logscore(z$logf[vi,,drop=FALSE],w[vi,,drop=FALSE])),0.)
    selected<-which.max(scores)
    dl<-lapply(ww,function(w)evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],
                z$nu[ti,,drop=FALSE],w[ti,,drop=FALSE]))
    delta_tables[[ticker]]<-do.call(rbind,lapply(seq_along(deltas),function(i) {
      e<-dl[[i]];data.frame(ticker=ticker,delta=deltas[i],validation_LPS=scores[i],
        selected_on_validation=i==selected,NLL=mean(e$NLL),QLIKE=mean(e$QLIKE),
        FZ0_025=mean(e$FZ_0p025))
    }))
    losses<-readRDS(file.path(outdir,"exploratory_full_pool",paste0(ticker,".rds")))$losses
    repair_tables[[ticker]]<-cbind(ticker=ticker,compare_losses(losses,"EWA_six",
      "EWA_six_old_transform",metrics=c("NLL","QLIKE","FZ_0p025"),B=999))
    losses$DMA_six_selected<-dl[[selected]]
    sg<-filter_compact(z$logf,simplex_lattice(6,5),0)$weights
    losses$Static_grid_six<-evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],
             z$nu[ti,,drop=FALSE],sg[ti,,drop=FALSE])
    for(m in c("Pool_filter_six","Pool_filter_six_censored","Static_grid_six","DMA_six_selected")) {
      ll<-losses[[m]]
      for(alpha in c(.01,.025,.05)) {
        s<-gsub("\\.","p",as.character(alpha));q<-ll[[paste0("VaR_",s)]];e<-ll[[paste0("ES_",s)]]
        calibration_tables[[length(calibration_tables)+1L]]<-cbind(ticker=ticker,method=m,alpha=alpha,
          coverage(z$y[ti],q,alpha),tail_calibration(z$y[ti],q,e,alpha))
      }
    }
    for(m in c("Pool_filter_six","Static_grid_six","Static_six")) {
      cc<-compare_losses(losses,m,"DMA_six_selected",metrics=c("NLL","FZ_0p025","QLIKE"),B=999)
      selected_comparisons[[length(selected_comparisons)+1L]]<-cbind(ticker=ticker,cc)
    }
    # Every unique comparison in the main and exploratory tables receives a
    # secondary QLIKE comparison. This is a separately disclosed test family.
    a<-read.csv(file.path(outdir,ticker,"comparisons.csv"))
    b<-read.csv(file.path(outdir,"exploratory_full_pool","comparisons.csv"))
    b<-b[b$ticker==ticker,]
    pairs<-unique(rbind(a[,c("candidate","comparator")],b[,c("candidate","comparator")]))
    for(i in seq_len(nrow(pairs))) {
      cc<-compare_losses(losses,pairs$candidate[i],pairs$comparator[i],metrics="QLIKE",B=999)
      qlike[[length(qlike)+1L]]<-cbind(ticker=ticker,cc)
    }
    # One declared mechanism contrast, all assets and three losses. Show both
    # lag and block sensitivity; do not select the most favorable specification.
    for(metric in c("NLL","QLIKE","FZ_0p025")) {
      d<-losses$Static_grid_six[[metric]]-losses$DMA_six_selected[[metric]]
      for(L in c(5,10,20,40)) {
        nw<-nw_difference(d,lag=L);ci<-block_interval(d,B=999,block=L)
        block_tables[[length(block_tables)+1L]]<-cbind(ticker=ticker,metric=metric,
                         block=L,nw,lo=ci[1],hi=ci[2])
      }
    }
    saveRDS(list(selected_delta=deltas[selected],selected_losses=dl[[selected]],
                static_grid_losses=losses$Static_grid_six,static_grid_weights=sg),
            file.path(dest,paste0(ticker,".rds")))
  }
  q<-do.call(rbind,qlike);q$p_holm_secondary_QLIKE<-p.adjust(q$p_two,"holm")
  c<-do.call(rbind,selected_comparisons);c$p_holm_selected_delta_family<-p.adjust(c$p_two,"holm")
  b<-do.call(rbind,block_tables);b$p_holm_robustness_family<-p.adjust(b$p_two,"holm")
  r<-do.call(rbind,repair_tables);r$p_holm_repair_family<-p.adjust(r$p_two,"holm")
  write.csv(do.call(rbind,delta_tables),file.path(dest,"six_pool_delta_sensitivity.csv"),row.names=FALSE)
  write.csv(q,file.path(dest,"secondary_QLIKE_comparisons.csv"),row.names=FALSE)
  write.csv(c,file.path(dest,"selected_delta_comparisons.csv"),row.names=FALSE)
  write.csv(b,file.path(dest,"block_and_lag_sensitivity.csv"),row.names=FALSE)
  write.csv(do.call(rbind,calibration_tables),file.path(dest,"extension_calibration.csv"),row.names=FALSE)
  write.csv(r,file.path(dest,"transformation_repair_comparisons.csv"),row.names=FALSE)
  invisible(dest)
}
