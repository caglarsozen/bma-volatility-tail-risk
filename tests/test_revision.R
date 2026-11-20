# Audit new outputs against raw returns, admissible training likelihoods and
# a direct EGARCH recursion. Mutation checks detect accidental target use.
run_revision_audit <- function(outdir = "results/diagnostics/verification",
                               tickers = default_config()$tickers) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cfg <- default_config()
  base <- "results/main_sample"
  eg <- "results/diagnostics/egarch_full"
  checks <- list()
  check <- function(name, error, tolerance = 0) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, maximum_error = error, tolerance = tolerance,
      pass = is.finite(error) && error <= tolerance)
  }
  for (ticker in tickers) {
    x <- read_prices(file.path("data/raw", paste0(ticker, ".csv")), ticker, cfg)
    b <- read.csv(file.path(base, ticker, "forecast_bank.csv"))
    for (nm in c("origin_date", "target_date", "fit_origin", "fit_target_end"))
      b[[nm]] <- as.Date(b[[nm]])
    check_bank(b, x, cfg)
    check(paste(ticker, "main bank matches frozen prices and information dates"), 0)

    a <- readRDS(file.path(eg, paste0(ticker, "_complete.rds")))
    z <- a$bank
    i <- match(z$target_date, x$date)
    first <- which(x$date >= cfg$validation_start)[1L]
    check(paste(ticker, "EGARCH alignment and completeness"), as.integer(
      !identical(i, first:nrow(x)) || any(z$origin_date != x$date[i-1L]) ||
      any(z$fit_target_end > z$origin_date) || any(!is.finite(z$sigma2)) ||
      any(z$sigma2 <= 0) || any(z$df <= 2)))
    check(paste(ticker, "EGARCH target values"), max(abs(z$ret - x$ret[i])), 1e-10)
    by_block <- split(a$diagnostics, a$diagnostics$block)
    gaps <- vapply(by_block, function(d) {
      if (sum(d$selected) != 1L || !d$valid[d$selected]) return(Inf)
      max(d$training_LL[d$valid]) - d$training_LL[d$selected]
    }, 0.)
    check(paste(ticker, "training maximum among admissible candidates"), max(abs(gaps)), 1e-10)
    blocks <- unique(as.integer(round(seq(1, length(by_block), length.out = 3))))
    mutation_error <- recursion_error <- 0
    for (block in blocks) {
      one <- readRDS(file.path(eg, paste0(ticker, "_block_", block, ".rds")))
      j <- match(one$rows$target_date, x$date)
      o <- one$n_train
      r <- x$ret[seq_len(tail(j, 1))]
      sp <- egarch_spec()
      rugarch::setfixed(sp) <- as.list(one$coefficients)
      v <- as.numeric(rugarch::sigma(rugarch::ugarchfilter(sp, r, n.old = o)))^2
      co <- one$coefficients
      nu <- co[["shape"]]
      abs_mean <- sqrt(nu - 2) * exp(lgamma((nu - 1)/2) - lgamma(nu/2)) / sqrt(pi)
      lagged_standardized <- r[j-1L] / sqrt(v[j-1L])
      implied <- exp(co[["omega"]] + co[["alpha1"]] * lagged_standardized +
        co[["gamma1"]] * (abs(lagged_standardized) - abs_mean) +
        co[["beta1"]] * log(v[j-1L]))
      recursion_error <- max(recursion_error, abs(v[j] - implied) / pmax(1, v[j]),
        abs(v[j] - one$rows$sigma2) / pmax(1, v[j]))
      pivot <- j[ceiling(length(j)/2)]
      r[pivot:length(r)] <- r[pivot:length(r)] + 11
      vm <- as.numeric(rugarch::sigma(rugarch::ugarchfilter(sp, r, n.old = o)))^2
      mutation_error <- max(mutation_error, abs(vm[j[j <= pivot]] - v[j[j <= pivot]]))
    }
    check(paste(ticker, "EGARCH direct one-step recursion"), recursion_error, 1e-7)
    check(paste(ticker, "EGARCH future-return mutation"), mutation_error, 1e-10)

    fr <- readRDS(file.path("results/diagnostics/refit5", ticker, "complete.rds"))$bank
    dates <- x$date[first:nrow(x)]
    expected <- first - 1L + 5L * ((seq_along(dates) - 1L) %/% 5L)
    for (m in unique(fr$model)) {
      rr <- fr[fr$model == m, ]
      rr <- rr[order(rr$target_date), ]
      check(paste(ticker, m, "five-observation training dates"), as.integer(
        !identical(rr$target_date, dates) || any(rr$fit_origin != x$date[expected]) ||
        any(rr$fit_target_end > rr$origin_date) || any(rr$origin_date >= rr$target_date)))
    }
  }
  ext <- read.csv("results/diagnostics/realized_benchmark/spyreal_official_scaled.csv")
  assert(requireNamespace("xts", quietly = TRUE), "xts required")
  data <- new.env()
  utils::data("spyreal", package = "rugarch", envir = data)
  target <- data$spyreal
  check("realized-kernel and return scaling from package data",
    max(abs(ext$r - 100 * as.numeric(target[, "SPY_OC"])),
        abs(ext$rv - 100 * as.numeric(target[, "SPY_RK"]))), 1e-10)
  dates <- as.Date(zoo::index(target))
  check("realized-kernel dates from package data",
    max(abs(as.numeric(as.Date(ext$date)) - as.numeric(dates))))
  result <- do.call(rbind, checks)
  write.csv(result, file.path(outdir, "revision_checks.csv"), row.names = FALSE)
  assert(all(result$pass), paste("Revision checks failed:", paste(result$check[!result$pass], collapse = ", ")))
  message(nrow(result), " new R checks passed")
  invisible(result)
}
