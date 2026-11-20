# Usage: Rscript 00_RUN.R audit | smoke | simulation | empirical | full
# Source this file to load functions only. No downloads, package installation or
# experiment is started by source(). Run from this project's root directory.
for(f in sort(list.files("R",pattern="\\.R$",full.names=TRUE)))source(f)
source("tests/test_core.R")
main <- function(args=commandArgs(trailingOnly=TRUE)) {
  mode<-if(length(args))args[1] else "audit"
  if(mode=="audit") return(run_tests("results/audit"))
  if(mode=="mechanism") return(run_log_proxy_mechanism())
  if(mode=="smoke")return(run_simulation("results/smoke",reps=2,n=800,validation_n=200))
  if(mode=="simulation") {
    reps<-as.integer(Sys.getenv("BMA_REPS","60"))
    return(run_simulation("results/simulation",reps=reps))
  }
  if(mode %in% c("empirical","full")) {
    if(mode=="full") {
      assert(!nchar(Sys.getenv("BMA_TICKERS","")),"Full mode requires all eight assets")
      run_tests("results/audit")
    }
    warning_log<-file.path("results",paste0("warnings_",format(Sys.time(),"%Y%m%d_%H%M%S"),".txt"))
    dir.create("results",showWarnings=FALSE)
    out<-withCallingHandlers(run_empirical(),warning=function(w) {
      cat(format(Sys.time()),conditionMessage(w),"\n",file=warning_log,append=TRUE)
      invokeRestart("muffleWarning")
    })
    if(mode=="full") {
      run_extra_audits(out);run_full_pool_extension(out);run_mechanism_audit(out)
      run_reviewer_robustness(out)
      run_log_proxy_mechanism()
    }
    return(out)
  }
  stop("Unknown mode: ",mode)
}
if(sys.nframe()==0L)main()
