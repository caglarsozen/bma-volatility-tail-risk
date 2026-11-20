# Forecasting volatility and tail risk under model uncertainty

R replication code for a study of variance scaling, complementary forecasts
and dynamic model averaging. The analysis compares equal-weight, fitted static and
adaptive pools of predictive densities, using volatility and joint tail-risk
losses on eight daily return series.

## Quick start

Run these commands from the repository root.

~~~sh
Rscript INSTALL_PACKAGES.R
Rscript 00_RUN.R audit
Rscript 00_RUN.R smoke
~~~

The package installer adds missing dependencies. The audit checks numerical
properties of the density and risk calculations. The smoke test runs a short
simulation, writing to `results/smoke/`. Both checks can run without market inputs.
The study used R 4.5.0 and `rugarch` 1.5-4.

## Market inputs

Place one comma-separated value (CSV) file per ticker in `data/raw/`. Files
are named `SPY.csv`, `QQQ.csv`, `IWM.csv`, `GLD.csv`, `TLT.csv`, `USO.csv`,
`BTC-USD.csv` and `ETH-USD.csv`. Each file has the following columns.

~~~text
date,open,high,low,close,adjusted,volume
~~~

Use `YYYY-MM-DD` for dates, the instrument's own trading calendar and the
adjusted closing price in `adjusted`. The code selects observations from
1 January 2018 through 29 June 2026. The validation period begins in 2021
and the evaluation period begins in 2022. The historical realized-kernel
case uses the `spyreal` example distributed with `rugarch`.

`Rscript DOWNLOAD_DATA.R` retrieves a new Yahoo Finance data vintage where
provider access permits. The script preserves existing snapshots. Use the
frozen snapshots from the study archive for the recorded analysis. Fresh
downloads and optimizer versions can change fitted forecasts.

## Full analysis

~~~sh
Rscript RUN_BMA_REVISION.R full
~~~

A full run creates a separate directory under `runs/`. It fits the
volatility models, evaluates the predictive pools and tail forecasts, and runs
the simulations. Each run writes the main forecasts to `results/main_sample/`
and the diagnostic analyses to `results/diagnostics/`. The completion marker
is `results/REVISION_COMPLETE.txt`. Resume an interrupted run in the same R and
package environment using its actual directory path.

~~~sh
Rscript RUN_BMA_REVISION.R resume runs/run_example
~~~

In RStudio, open the repository root and run

~~~r
source("RUN_BMA_REVISION.R")
bma_revision("full")
~~~

## Verify existing results

Run `Rscript RUN_BMA_REVISION.R verify` from a completed run directory
containing the price snapshots and forecast outputs. For a study archive,
use the verification scripts bundled with that archive.

The independent Python audits in `tests/` use NumPy, pandas and SciPy.
`AUDIT_EXPORT.R` exports per-date losses and outcomes for
`tests/independent_audit.py` when run alongside the saved forecasts.

| Path | Purpose |
|---|---|
| `R/` | Density pools, volatility models, simulation and empirical analyses |
| `tests/` | Numerical and forecast chronology checks |
| `DOWNLOAD_DATA.R` | Retrieve a separate market-data vintage |
| `RUN_BMA_REVISION.R` | Run, resume or verify the full study |

For a Turkish quick-start guide, see [R çalıştırma yönergesi](docs/RUN_TR.md).

## License

The project code is released under the [MIT License](LICENSE). Market inputs
obtained from data providers and examples distributed with external packages
retain their respective terms.
