#!/usr/bin/env Rscript
# Independent populations can use separate R processes; each output directory
# has exactly one writer. The final serial call aggregates validated caches.
kit <- normalizePath('benchmark/time_series_jasinska2020')
args <- commandArgs(TRUE)
workers <- as.integer(if(length(args))args[1] else Sys.getenv('SLURM_CPUS_PER_TASK','2'))
stopifnot(length(workers)==1,!is.na(workers),workers>=1)
wells <- unique(read.delim(file.path(kit,'samples.tsv'))$well)
workers <- min(workers,length(wells))
logs <- file.path(kit,'generated/clustering/logs');dir.create(logs,recursive=TRUE,showWarnings=FALSE)
script <- file.path(kit,'cluster_time_series.R')
# Verify the separately installed optimized package before launching workers.
source(file.path(kit,'load_release.R'));invisible(load_release(kit))
status <- parallel::mclapply(wells,function(w) {
  log <- file.path(logs,paste0(w,'.log'))
  system2(file.path(R.home('bin'),'Rscript'),c(shQuote(script),w),stdout=log,stderr=log)
},mc.cores=workers)
ok <- vapply(status,function(s)is.numeric(s) && length(s)==1 && !is.na(s) && s==0,logical(1))
if(!all(ok))stop('Clustering failed for: ',paste(wells[!ok],collapse=', '),'. See ',logs)
stopifnot(system2(file.path(R.home('bin'),'Rscript'),shQuote(script))==0)
