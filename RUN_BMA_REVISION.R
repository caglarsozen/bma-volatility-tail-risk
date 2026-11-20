# Rscript RUN_BMA_REVISION.R full
# Rscript RUN_BMA_REVISION.R resume <run-directory>
# Rscript RUN_BMA_REVISION.R verify
# From RStudio at the project root: source("RUN_BMA_REVISION.R"); bma_revision("full")

revision_manifest <- function(directory) {
  rel <- c("00_RUN.R", "RUN_BMA_REVISION.R",
           sort(list.files(file.path(directory, "R"), pattern = "\\.R$")),
           sort(list.files(file.path(directory, "data", "raw"), pattern = "\\.csv$")),
           sort(list.files(file.path(directory, "tests"), pattern = "\\.R$")))
  prefix <- c("", "", rep("R", length(list.files(file.path(directory, "R"), pattern = "\\.R$"))),
              rep("data/raw", length(list.files(file.path(directory, "data", "raw"), pattern = "\\.csv$"))),
              rep("tests", length(list.files(file.path(directory, "tests"), pattern = "\\.R$"))))
  files <- ifelse(prefix == "", rel, paste(prefix, rel, sep = "/"))
  hashes <- unname(tools::md5sum(file.path(directory, files)))
  stopifnot(!anyNA(hashes))
  hashes <- setNames(hashes, files)
  packages <- c("rugarch", "Rsolnp", "xts", "nloptr", "numDeriv")
  list(hashes = hashes, R = R.version.string, platform = R.version$platform,
       packages = vapply(packages, function(x) as.character(utils::packageVersion(x)), ""))
}

bma_revision <- function(mode = "full", run_directory = NULL) {
  mode <- match.arg(mode, c("full", "resume", "verify"))
  root <- getwd()
  on.exit(setwd(root), add = TRUE)
  if (!file.exists("00_RUN.R")) stop("Run from the repository root containing 00_RUN.R.")
  if (mode != "resume") {
    source("00_RUN.R")
    tickers <- default_config()$tickers
    snapshots <- file.path("data", "raw", paste0(tickers, ".csv"))
    missing <- snapshots[!file.exists(snapshots)]
    if (length(missing)) stop(
      "Missing price snapshots: ", paste(missing, collapse = ", "),
      ". See the market-data setup in README.md.", call. = FALSE)
  }
  for (pkg in c("rugarch", "xts", "nloptr", "numDeriv")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install required package: ", pkg)
  }
  if (mode == "verify") {
    baseline <- file.path("results", "main_sample")
    missing <- !dir.exists(baseline) || !dir.exists("results/diagnostics/egarch_full") ||
      !dir.exists("results/diagnostics/refit5") || !dir.exists("results/diagnostics/realized_benchmark")
    if (missing) stop(
      "Archived forecast outputs are required by verify mode. See README.md.",
      call. = FALSE)
    source("tests/test_revision.R")
    run_tests("results/diagnostics/verification_core")
    run_revision_audit()
    return(invisible(root))
  }

  if (mode == "full") {
    if (is.null(run_directory)) run_directory <- file.path(
      "runs", paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
    if (file.exists(run_directory)) stop("Full mode requires a new output directory.")
    dir.create(run_directory, recursive = TRUE)
    dir.create(file.path(run_directory, "data"))
    stopifnot(all(file.copy(c("00_RUN.R", "RUN_BMA_REVISION.R", "R", "tests"),
      run_directory, recursive = TRUE)))
    stopifnot(file.copy("data/raw", file.path(run_directory, "data"), recursive = TRUE))
    saveRDS(revision_manifest(run_directory), file.path(run_directory, "revision_inputs.rds"))
  } else {
    if (is.null(run_directory)) stop("Give the output directory of the interrupted run.")
    receipt <- file.path(run_directory, "revision_inputs.rds")
    if (!file.exists(receipt)) stop("The run directory has no revision input record.")
    recorded <- readRDS(receipt)
    if (!identical(recorded, revision_manifest(run_directory)))
      stop("The input bytes or software versions differ from the interrupted run.")
    current_runner <- unname(tools::md5sum("RUN_BMA_REVISION.R"))
    if (!identical(current_runner, recorded$hashes[["RUN_BMA_REVISION.R"]])) {
      saved_directory <- normalizePath(run_directory, mustWork = TRUE)
      setwd(saved_directory)
      saved_runner <- new.env()
      sys.source("RUN_BMA_REVISION.R", envir = saved_runner)
      return(saved_runner$bma_revision("resume", saved_directory))
    }
  }
  run_directory <- normalizePath(run_directory, mustWork = TRUE)
  setwd(run_directory)
  source("00_RUN.R")
  source("tests/test_revision.R")
  dir.create("results", showWarnings = FALSE)
  cfg <- default_config()
  set.seed(cfg$seed)
  sys <- Sys.getenv(c("BMA_RUN_ID", "BMA_TICKERS", "BMA_BOOTSTRAP"), unset = NA_character_)
  on.exit({
    if (any(!is.na(sys))) do.call(Sys.setenv, as.list(sys[!is.na(sys)]))
    Sys.unsetenv(names(sys)[is.na(sys)])
  }, add = TRUE)
  Sys.setenv(BMA_RUN_ID = "main_sample", BMA_TICKERS = "", BMA_BOOTSTRAP = "999")

  run_tests("results/audit")
  run_egarch_multistart_bank()
  run_egarch_curvature_audit()
  refine_egarch_curvature()
  fresh_volatility <- function(x, config) {
    pair <- garch_bank(x, config, models = c(GARCH = "sGARCH", GJR = "gjrGARCH"))
    repaired <- readRDS(file.path("results/diagnostics/egarch_full", paste0(x$ticker[1], "_complete.rds")))$bank
    rbind(pair, repaired)
  }
  # The fresh main bank starts with repaired EGARCH. The supplied archive
  # retains the original Windows forecasts for the numerical before/after comparison.
  daily <- run_empirical(cfg, volatility_builder = fresh_volatility, cache_tag = "repaired-egarch")
  run_extra_audits(daily)
  run_full_pool_extension(daily)
  run_mechanism_audit(daily)
  run_reviewer_robustness(daily)
  run_log_proxy_mechanism()
  run_proxy_granularity()
  run_realized_measure_benchmark()
  run_calibration_and_floor()
  run_scale_weight_attribution()
  run_repaired_pool_evaluation(record_legacy_ranges = FALSE)
  run_refit_sensitivity()
  collect_refit_sensitivity()
  if (!file.exists("results/simulation/simulation_comparisons.csv"))
    run_simulation("results/simulation", reps = 60)
  run_revision_audit()
  writeLines("Completed eight-asset revision from frozen prices", "results/REVISION_COMPLETE.txt")
  message("Completed: ", run_directory)
  invisible(run_directory)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  path <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  supplied <- if (length(args) >= 2L) args[[2]] else NULL
  if (identical(args[1L], "resume") && !is.null(supplied))
    supplied <- normalizePath(supplied, mustWork = TRUE)
  if (length(path) == 1L) setwd(dirname(normalizePath(sub("^--file=", "", path))))
  bma_revision(if (length(args)) args[[1]] else "full", supplied)
}
