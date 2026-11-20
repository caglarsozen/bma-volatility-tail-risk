source("00_RUN.R")
outdir<-"results/main_sample"
dir.create("results/independent_audit",showWarnings=FALSE)
for(ticker in default_config()$tickers) {
  ll<-readRDS(file.path(outdir,ticker,"daily_losses.rds"))
  for(m in c("Pool_filter","Static_pair","DMA_six","EWA_six","EWA_six_old_transform")) {
    z<-ll[[m]]
    write.csv(z,file.path("results/independent_audit",paste0(ticker,"_",m,".csv")),row.names=FALSE)
  }
}
run_extra_audits(outdir)
