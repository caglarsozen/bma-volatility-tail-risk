# Explicit, optional package installation. Never called by 00_RUN.R.
packages<-c("rugarch","xts","jsonlite","nloptr","numDeriv")
missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos="https://cloud.r-project.org")
