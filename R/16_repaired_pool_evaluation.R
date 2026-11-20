# Replace only the EGARCH rows of the user's recorded forecast bank. All six
# original forecasts, raw targets and dates are held fixed for this comparison.
run_repaired_pool_evaluation <- function(outdir="results/diagnostics/repaired_pool",
                                         record_legacy_ranges=TRUE,
                                         tickers=default_config()$tickers) {
  cfg<-default_config();dir.create(outdir,recursive=TRUE,showWarnings=FALSE);base<-"results/main_sample"
  perf<-comp<-selection<-ranges<-list()
  for(ticker in tickers) {
    dest<-file.path(outdir,paste0(ticker,"_complete.rds"))
    if(file.exists(dest))a<-readRDS(dest) else {
      b<-read.csv(file.path(base,ticker,"forecast_bank.csv"));for(nm in c("origin_date","target_date","fit_origin","fit_target_end"))b[[nm]]<-as.Date(b[[nm]])
      eb<-readRDS(file.path("results/diagnostics/egarch_full",paste0(ticker,"_complete.rds")))$bank
      old<-b[b$model=="EGARCH",];old<-old[order(old$target_date),];eb<-eb[order(eb$target_date),]
      assert(identical(old$target_date,eb$target_date)&&max(abs(old$ret-eb$ret))<1e-10,"EGARCH replacement target mismatch")
      a<-evaluate_revision_pool(rbind(b[b$model!="EGARCH",],eb),cfg,c("GARCH","GJR","EGARCH","HAR","Range_HAR","HAR66","Log_AR1"))
      z<-bank_matrices(eb,"EGARCH");ti<-which(z$dates>=cfg$test_start)
      a$losses$EGARCH_repaired<-evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],matrix(1,length(ti),1))
      a$performance<-rbind(a$performance,data.frame(method="EGARCH_repaired",n=length(ti),NLL=mean(a$losses$EGARCH_repaired$NLL),QLIKE=mean(a$losses$EGARCH_repaired$QLIKE),FZ0_025=mean(a$losses$EGARCH_repaired$FZ_0p025)))
      orig<-readRDS(file.path(base,ticker,"daily_losses.rds"));s<-readRDS(file.path(base,"reviewer_robustness",paste0(ticker,".rds")))
      a$losses$DMA_six_selected<-s$selected_losses;a$losses$Static_six<-orig$Static_six;a$losses$Static_grid_six<-s$static_grid_losses
      a$comparisons<-do.call(rbind,lapply(c("DMA_six_selected","Static_six","Static_grid_six"),function(m)
        compare_losses(a$losses,m,c("Static_pair","EGARCH_repaired","Static","DMA_selected","Filter"),B=999)))
      a$ranges<-data.frame(ticker=ticker,max_original_variance=max(old$sigma2),max_repaired_variance=max(eb$sigma2),
        original_over_1e6=sum(old$sigma2>1e6),repaired_over_1e6=sum(eb$sigma2>1e6),min_repaired_df=min(eb$df))
      saveRDS(a,dest)
    }
    perf[[ticker]] <- cbind(ticker = ticker, a$performance)
    cc <- a$comparisons
    cc$maximum_daily_loss_difference <- vapply(seq_len(nrow(cc)), function(j) {
      metric <- cc$metric[j]
      max(abs(a$losses[[cc$candidate[j]]][[metric]] -
              a$losses[[cc$comparator[j]]][[metric]]))
    }, 0.)
    # Keep the raw tests. Machine-scale differences between two independently
    # optimized weights can generate a large standardized statistic.
    cc$numerically_tied <- cc$maximum_daily_loss_difference <= 1e-8
    comp[[ticker]] <- cbind(ticker = ticker, cc)
    selection[[ticker]]<-cbind(ticker=ticker,a$selection);ranges[[ticker]]<-a$ranges
    message("Repaired seven-pool evaluation: ",ticker)
  }
  if(!identical(tickers,cfg$tickers))return(invisible(outdir))
  c<-do.call(rbind,comp);c$p_holm_repaired_pool_family<-p.adjust(c$p_two,"holm",n=nrow(c))
  write.csv(do.call(rbind,perf),file.path(outdir,"performance.csv"),row.names=FALSE)
  write.csv(c,file.path(outdir,"comparisons.csv"),row.names=FALSE)
  write.csv(do.call(rbind,selection),file.path(outdir,"selection.csv"),row.names=FALSE)
  if(record_legacy_ranges)write.csv(do.call(rbind,ranges),file.path(outdir,"variance_ranges.csv"),row.names=FALSE)
  invisible(outdir)
}
