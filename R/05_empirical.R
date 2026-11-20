bank_matrices <- function(bank,models) {
  dates<-sort(unique(bank$target_date));sub<-bank[bank$model==models[1],];sub<-sub[order(sub$target_date),]
  assert(identical(sub$target_date,dates),"Reference bank incomplete")
  v<-nu<-matrix(NA_real_,length(dates),length(models),dimnames=list(NULL,models))
  for(j in seq_along(models)) {
    x<-bank[bank$model==models[j],];x<-x[order(x$target_date),]
    assert(identical(x$target_date,dates)&&max(abs(x$ret-sub$ret))<1e-12,"Unequal model targets")
    v[,j]<-x$sigma2;nu[,j]<-x$df
  }
  list(dates=dates,origins=sub$origin_date,y=sub$ret,v=v,nu=nu,logf=t_logpdf(sub$ret,v,nu))
}
evaluate_bank <- function(bank,config,outdir,B=999) {
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  z<-bank_matrices(bank,c("GARCH","GJR"));n<-length(z$y)
  vi<-which(z$dates<config$test_start);ti<-which(z$dates>=config$test_start)
  cut<-mix_risk(matrix(.5,n,2),z$v,z$nu,.10)$q
  cscore<-censored_logf(z$y,z$logf,z$v,z$nu,cut)
  methods<-build_methods(z$logf,vi,cscore)
  losses<-list();weights<-list()
  for(m in names(methods$weights)) {
    w<-methods$weights[[m]];e<-evaluate_mixture(z$y[ti],z$v[ti,,drop=FALSE],z$nu[ti,,drop=FALSE],w[ti,,drop=FALSE])
    e$date<-z$dates[ti];e$ret<-z$y[ti];losses[[m]]<-e
    weights[[m]]<-data.frame(date=z$dates,method=m,GARCH=w[,1],GJR=w[,2])
  }
  # Strong single models and corrected six-/seven-model pool controls.
  sets<-list(six=c("GARCH","GJR","HAR","Range_HAR","HAR66","Log_AR1"),
             seven=c("GARCH","GJR","EGARCH","HAR","Range_HAR","HAR66","Log_AR1"))
  for(pool in names(sets)) {
    a<-bank_matrices(bank,sets[[pool]]);nm<-length(sets[[pool]])
    sw<-static_pool(a$logf[vi,,drop=FALSE]);ew<-matrix(1/nm,n,nm)
    ww<-list(EWA=ew,Static=matrix(rep(sw,each=n),n,nm),DMA=dma_weights(a$logf,.97))
    for(m in names(ww)) {
      e<-evaluate_mixture(a$y[ti],a$v[ti,,drop=FALSE],a$nu[ti,,drop=FALSE],ww[[m]][ti,,drop=FALSE])
      e$date<-a$dates[ti];e$ret<-a$y[ti];losses[[paste(m,pool,sep="_")]]<-e
    }
    if(pool=="six") {
      legacy<-a
      for(j in 3:6) {
        b<-bank[bank$model==sets[[pool]][j],];b<-b[order(b$target_date),]
        legacy$v[,j]<-b$sigma2_lognormal
      }
      e<-evaluate_mixture(a$y[ti],legacy$v[ti,,drop=FALSE],a$nu[ti,,drop=FALSE],ew[ti,,drop=FALSE])
      e$date<-a$dates[ti];e$ret<-a$y[ti];losses$EWA_six_old_transform<-e
    }
  }
  for(m in unique(bank$model)) {
    a<-bank_matrices(bank,m);e<-evaluate_mixture(a$y[ti],a$v[ti,,drop=FALSE],a$nu[ti,,drop=FALSE],matrix(1,length(ti),1))
    e$date<-a$dates[ti];e$ret<-a$y[ti];losses[[m]]<-e
  }
  perf<-do.call(rbind,lapply(names(losses),function(m) {
    l<-losses[[m]];data.frame(method=m,n=nrow(l),LPS=-mean(l$NLL),NLL=mean(l$NLL),
      QLIKE=mean(l$QLIKE),FZ0_025=mean(l$FZ_0p025))
  }))
  comparisons<-do.call(rbind,lapply(c("Pool_filter","Pool_filter_censored","DMA_six"),function(m)
    compare_losses(losses,m,c("EWA_pair","Static_pair","GARCH","GJR","DMA_selected","LDF2","Rolling_pool","Static_censored"),B=B)))
  # Explicit family: all three candidates x eight comparators x two metrics per asset.
  comparisons$p_holm_asset_family<-p.adjust(comparisons$p_two,"holm")
  bt<-list();pit<-list();k<-0L
  for(m in names(losses)) {
    l<-losses[[m]];pit[[m]]<-cbind(method=m,pit_diagnostics(l$PIT))
    for(alpha in c(.01,.025,.05)) {
      suffix<-gsub("\\.","p",as.character(alpha));q<-l[[paste0("VaR_",suffix)]];e<-l[[paste0("ES_",suffix)]]
      k<-k+1;bt[[k]]<-cbind(method=m,alpha=alpha,coverage(l$ret,q,alpha),tail_calibration(l$ret,q,e,alpha))
    }
  }
  # Exact, additive score accounting on the same observations, not causal effects.
  ladder<-c("EWA_six_old_transform","EWA_six","EWA_pair","Static_pair","Pool_filter")
  decomp<-do.call(rbind,lapply(c("NLL","FZ_0p025","QLIKE"),function(metric) {
    av<-vapply(ladder,function(m)mean(losses[[m]][[metric]]),0.)
    data.frame(metric=metric,step=c("HAR_transformation","pool_restriction","static_learning","dynamic_pool_learning"),
      from=head(ladder,-1),to=tail(ladder,-1),loss_reduction=head(av,-1)-tail(av,-1))
  }))
  # Tail-shape sensitivity: fixed mapping only, no re-estimation or test selection.
  sensitivity<-list();kk<-0
  for(df in c(5,8,12)) {
    a<-bank_matrices(bank,sets$six);a$nu[,3:6]<-df;a$logf<-t_logpdf(a$y,a$v,a$nu)
    for(method in c("EWA_six","DMA_six")) {
      w<-if(method=="EWA_six")matrix(1/6,n,6) else dma_weights(a$logf,.97)
      e<-evaluate_mixture(a$y[ti],a$v[ti,,drop=FALSE],a$nu[ti,,drop=FALSE],w[ti,,drop=FALSE],alphas=.025)
      kk<-kk+1;sensitivity[[kk]]<-data.frame(method=method,HAR_df=df,LPS=-mean(e$NLL),QLIKE=mean(e$QLIKE),FZ0_025=mean(e$FZ_0p025))
    }
  }
  write.csv(perf,file.path(outdir,"performance.csv"),row.names=FALSE)
  write.csv(comparisons,file.path(outdir,"comparisons.csv"),row.names=FALSE)
  write.csv(do.call(rbind,bt),file.path(outdir,"calibration.csv"),row.names=FALSE)
  write.csv(do.call(rbind,pit),file.path(outdir,"pit.csv"),row.names=FALSE)
  write.csv(do.call(rbind,weights),file.path(outdir,"pair_forecast_weights.csv"),row.names=FALSE)
  write.csv(methods$selection,file.path(outdir,"validation_selection.csv"),row.names=FALSE)
  write.csv(decomp,file.path(outdir,"gain_accounting.csv"),row.names=FALSE)
  write.csv(do.call(rbind,sensitivity),file.path(outdir,"df_sensitivity.csv"),row.names=FALSE)
  saveRDS(losses,file.path(outdir,"daily_losses.rds"))
  invisible(list(performance=perf,comparisons=comparisons))
}
run_empirical <- function(config=default_config(), volatility_builder=garch_bank, cache_tag="baseline") {
  set.seed(config$seed)
  chosen<-Sys.getenv("BMA_TICKERS","");if(nchar(chosen))config$tickers<-strsplit(chosen,",",fixed=TRUE)[[1]]
  runid<-Sys.getenv("BMA_RUN_ID",paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  outdir<-file.path("results",runid);dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  # Reuse allowed only for matching input bytes, exact configuration, and model source.
  cache<-file.path("data","forecast_cache");dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  perf<-comps<-quality<-list()
  for(ticker in config$tickers) {
    path<-file.path("data","raw",paste0(ticker,".csv"));assert(file.exists(path),paste("Missing raw snapshot",path))
    x<-read_prices(path,ticker,config)
    key<-list(config=config[setdiff(names(config),"tickers")],ticker=ticker,
      md5=unname(tools::md5sum(c(path,"R/03_models.R"))))
    if(cache_tag!="baseline")key$forecast_builder<-cache_tag
    dest<-file.path(cache,paste0(ticker,".rds"));cached<-if(file.exists(dest))readRDS(dest) else NULL
    if(!is.null(cached)&&identical(cached$key,key)) {bank<-cached$bank;message("Using checked forecast cache: ",ticker)} else {
      bank<-rbind(har_bank(x,config),volatility_builder(x,config));check_bank(bank,x,config)
      saveRDS(list(key=key,bank=bank),dest)
    }
    check_bank(bank,x,config)
    td<-file.path(outdir,ticker);dir.create(td,recursive=TRUE,showWarnings=FALSE)
    write.csv(bank,file.path(td,"forecast_bank.csv"),row.names=FALSE)
    z<-evaluate_bank(bank,config,td,B=as.integer(Sys.getenv("BMA_BOOTSTRAP","999")))
    perf[[ticker]]<-cbind(ticker=ticker,z$performance);comps[[ticker]]<-cbind(ticker=ticker,z$comparisons)
    quality[[ticker]]<-data.frame(ticker=ticker,first=x$date[1],last=tail(x$date,1),n=nrow(x),
      mean_return=mean(x$ret),max_abs_return=max(abs(x$ret)),zero_returns=sum(x$ret==0),
      log_floor_count=sum(x$rv<config$log_floor),calendar="native",validation_n=sum(x$date>=config$validation_start&x$date<config$test_start),
      evaluation_n=sum(x$date>=config$test_start),source_vintage="local_snapshot")
  }
  pp<-do.call(rbind,perf);cc<-do.call(rbind,comps)
  cc$p_holm_global<-p.adjust(cc$p_two,"holm")
  write.csv(pp,file.path(outdir,"all_performance.csv"),row.names=FALSE)
  write.csv(cc,file.path(outdir,"all_comparisons.csv"),row.names=FALSE)
  write.csv(do.call(rbind,quality),file.path(outdir,"data_audit.csv"),row.names=FALSE)
  saveRDS(config,file.path(outdir,"config.rds"))
  writeLines(capture.output(dput(config)),file.path(outdir,"config.txt"))
  writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
  writeLines("COMPLETED: all requested assets and seven models passed completeness checks",file.path(outdir,"STATUS.txt"))
  invisible(outdir)
}
