# Explicit download only. Existing snapshots are preserved; use a new directory
# or manually archive them before obtaining a new vintage.
source("R/01_density_pool.R");source("R/03_models.R")
download_snapshots <- function(config=default_config()) {
  assert(requireNamespace("jsonlite",quietly=TRUE),"Install jsonlite first")
  dir.create("data/raw",recursive=TRUE,showWarnings=FALSE)
  metadata<-list()
  for(ticker in config$tickers) {
    path<-file.path("data/raw",paste0(ticker,".csv"))
    assert(!file.exists(path),paste("Snapshot exists and will not be overwritten:",path))
    stamp<-function(d)as.numeric(as.POSIXct(d,tz="UTC"))
    url<-sprintf("https://query2.finance.yahoo.com/v8/finance/chart/%s?period1=%.0f&period2=%.0f&interval=1d",ticker,stamp(config$start),stamp(config$end_exclusive))
    jp<-file.path("data/raw",paste0(ticker,".json"))
    download.file(url,jp,quiet=TRUE,headers=c(`User-Agent`="Mozilla/5.0"))
    j<-jsonlite::fromJSON(jp,simplifyVector=FALSE)$chart$result[[1]]
    assert(!is.null(j),"Yahoo returned no data")
    num<-function(z)vapply(z,function(v)if(is.null(v))NA_real_ else as.numeric(v),0.)
    quote<-j$indicators$quote[[1]];tz<-j$meta$exchangeTimezoneName
    date<-as.Date(as.POSIXct(num(j$timestamp),origin="1970-01-01",tz=tz),tz=tz)
    d<-data.frame(date=date,open=num(quote$open),high=num(quote$high),low=num(quote$low),
      close=num(quote$close),adjusted=num(j$indicators$adjclose[[1]]$adjclose),volume=num(quote$volume))
    write.csv(d,path,row.names=FALSE)
    metadata[[ticker]]<-data.frame(ticker=ticker,downloaded_utc=format(Sys.time(),tz="UTC",usetz=TRUE),timezone=tz,
      source=url,raw_md5=unname(tools::md5sum(jp)),vintage="new_download_not_original_snapshot")
  }
  write.csv(do.call(rbind,metadata),"data/raw/download_manifest.csv",row.names=FALSE)
}
if(sys.nframe()==0L)download_snapshots()
