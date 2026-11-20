# Paired 2 x 2 experiment with identical dates, regression fits and t shapes
# within each panel. The only switches are the back-transformation and weights.
scale_weight_cells <- function(z, old_variances, evaluation, delta = .97) {
  assert(identical(dim(z$v), dim(old_variances)), "Variance panels differ")
  equal <- matrix(1 / ncol(z$v), nrow(z$v), ncol(z$v))
  losses <- list()
  for (transform in c("Smearing", "Normal")) {
    v <- if (transform == "Smearing") z$v else old_variances
    dma <- dma_weights(t_logpdf(z$y, v, z$nu), delta)
    for (mix in c("Equal", "DMA")) {
      w <- if (mix == "Equal") equal else dma
      losses[[paste(transform, mix, sep = "_")]] <- evaluate_mixture(
        z$y[evaluation], v[evaluation, , drop = FALSE],
        z$nu[evaluation, , drop = FALSE], w[evaluation, , drop = FALSE],
        alphas = .025
      )
    }
  }
  losses
}

scale_weight_interaction <- function(losses, metric) {
  (losses$Normal_Equal[[metric]] - losses$Normal_DMA[[metric]]) -
    (losses$Smearing_Equal[[metric]] - losses$Smearing_DMA[[metric]])
}

scale_weight_summary <- function(losses) {
  metrics <- c("NLL", "QLIKE", "FZ_0p025")
  do.call(rbind, lapply(metrics, function(metric) {
    a <- losses$Normal_Equal[[metric]]
    b <- losses$Normal_DMA[[metric]]
    c <- losses$Smearing_Equal[[metric]]
    d <- losses$Smearing_DMA[[metric]]
    # Positive means the old transformation amplifies the apparent DMA gain.
    gain_change <- scale_weight_interaction(losses, metric)
    nw <- nw_difference(gain_change)
    ci <- block_interval(gain_change, B = 999, block = 20)
    data.frame(metric = metric, old_equal = mean(a), old_DMA = mean(b),
      corrected_equal = mean(c), corrected_DMA = mean(d),
      old_advantage = mean(a - b), corrected_advantage = mean(c - d),
      nw, lo = ci[1], hi = ci[2])
  }))
}

run_scale_weight_attribution <- function(outdir = "results/diagnostics/scale_weight") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cfg <- default_config()
  base <- "results/main_sample"
  models <- c("GARCH", "GJR", "HAR", "Range_HAR", "HAR66", "Log_AR1")
  main <- list()
  for (ticker in cfg$tickers) {
    bank <- read.csv(file.path(base, ticker, "forecast_bank.csv"))
    bank$target_date <- as.Date(bank$target_date)
    z <- bank_matrices(bank, models)
    old_v <- z$v
    for (j in 3:6) {
      a <- bank[bank$model == models[j], ]
      old_v[, j] <- a$sigma2_lognormal[order(a$target_date)]
    }
    losses <- scale_weight_cells(z, old_v, z$dates >= cfg$test_start)
    main[[ticker]] <- cbind(ticker = ticker, scale_weight_summary(losses))
  }
  main <- do.call(rbind, main)
  main$p_holm_attribution_family <- p.adjust(main$p_two, method = "holm", n = nrow(main))
  write.csv(main, file.path(outdir, "main_factorial.csv"), row.names = FALSE)

  # Compare measurement constructions on the same SPY open-to-close dates.
  # This is an observed-proxy comparison, not a randomized sampling intervention.
  external <- "results/diagnostics/realized_benchmark"
  obj <- readRDS(file.path(external, "complete.rds"))
  x <- read.csv(file.path(external, "spyreal_official_scaled.csv"))
  x$date <- as.Date(x$date)
  x$ret <- x$r
  x$rv <- x$r^2
  x$rs <- x$rv
  x$ticker <- "SPY_external"
  extcfg <- cfg
  extcfg$validation_start <- as.Date("2004-01-01")
  extcfg$test_start <- as.Date("2005-01-01")
  extcfg$min_train <- 400L
  # Only the HAR(1,5,22) rows enter the squared-return design.
  b <- har_bank(x, extcfg)
  common <- obj$bank[obj$bank$model != "HAR_RK", ]
  har <- b[b$model == "HAR", names(common)]
  ext_rows <- ext_losses <- list()
  for (proxy in c("Realized kernel", "Squared open-to-close return")) {
    b <- if (proxy == "Realized kernel") obj$bank else rbind(common, har)
    model <- if (proxy == "Realized kernel") "HAR_RK" else "HAR"
    z <- bank_matrices(b, c("GARCH", "GJR", "RealGARCH", model))
    old <- z$v
    m <- b[b$model == model, ]
    m <- m[order(m$target_date), ]
    old[, 4] <- m$sigma2_lognormal
    eligible <- z$dates >= extcfg$test_start
    assert(identical(z$dates[eligible], obj$dates), "External dates differ")
    losses <- scale_weight_cells(z, old, eligible)
    ext_losses[[proxy]] <- losses
    ext_rows[[proxy]] <- cbind(proxy = proxy,
      median_normal_smear_ratio = median(old[, 4] / z$v[, 4]),
      scale_weight_summary(losses))
    write.csv(m, file.path(outdir, if (model == "HAR_RK") "external_HAR_RK.csv" else "external_HAR_R2.csv"), row.names = FALSE)
  }
  ext <- do.call(rbind, ext_rows)
  ext$p_holm_external_attribution <- p.adjust(ext$p_two, "holm", n = nrow(ext))
  write.csv(ext, file.path(outdir, "external_factorial.csv"), row.names = FALSE)
  # A direct paired contrast, rather than comparing two significance indicators.
  contrast <- do.call(rbind, lapply(c("NLL", "QLIKE", "FZ_0p025"), function(metric) {
    d <- scale_weight_interaction(ext_losses[["Squared open-to-close return"]], metric) -
      scale_weight_interaction(ext_losses[["Realized kernel"]], metric)
    ci <- block_interval(d, B = 999, block = 20)
    data.frame(metric = metric, nw_difference(d), lo = ci[1], hi = ci[2])
  }))
  contrast$p_holm_proxy_contrast <- p.adjust(contrast$p_two, "holm")
  write.csv(contrast, file.path(outdir, "external_paired_proxy_contrast.csv"), row.names = FALSE)
  invisible(outdir)
}
