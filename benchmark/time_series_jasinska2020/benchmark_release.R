#!/usr/bin/env Rscript
# Re-time existing development results with an optimized binary and verify every
# membership and centroid count. A matching release run can reuse its own timing.
kit <- normalizePath('benchmark/time_series_jasinska2020')
source(file.path(kit,'load_release.R'));build <- load_release(kit)
source(file.path(kit,'provenance_helpers.R'))
sha <- function(f)digest::digest(file=f,algo='sha256')
args <- commandArgs(TRUE);workers <- if(length(args))as.integer(args[1]) else 3L
stopifnot(length(workers)==1,!is.na(workers),workers>=1)
wells <- unique(read.delim(file.path(kit,'samples.tsv'))$well)
workers <- min(workers,length(wells))
dest <- file.path(kit,'generated/release/checks');dir.create(dest,showWarnings=FALSE)
check <- function(well) {
  cdir <- file.path(kit,'generated/clustering',well)
  deadline <- Sys.time()+7200
  while(!file.exists(file.path(cdir,'completed.rds'))) {
    if(Sys.time()>deadline)stop('Missing complete clustering: ',well)
    Sys.sleep(10)
  }
  files <- file.path(cdir,c('input.csv','members.csv','centroids.csv','completed.rds'))
  hashes <- setNames(vapply(files,sha,character(1)),basename(files))
  signature <- release_signature(hashes,build)
  cache <- file.path(dest,paste0(well,'.rds'))
  if(file.exists(cache)) {
    old <- readRDS(cache)
    if(identical(old$signature,signature))return(old)
  }
  previous <- readRDS(file.path(cdir,'completed.rds'))
  if(identical(previous$binary_sha256,build$binary_sha256)) {
    tm <- previous$timing
    mode <- 'Original clustering used this optimized binary'
  } else {
    input <- readr::read_csv(files[1],show_col_types=FALSE)
    elapsed <- system.time(calls <- do.call(barbac::super_cluster2,
      c(list(input),release_settings(),list(verbose=TRUE))))
    members <- readr::read_csv(files[2],show_col_types=FALSE)
    centroids <- readr::read_csv(files[3],show_col_types=FALSE)
    barcode <- unlist(calls$all_barcodes,use.names=FALSE)
    centroid <- rep(calls$central_barcode,lengths(calls$all_barcodes))
    stopifnot(!anyDuplicated(barcode),setequal(barcode,members$barcode),
      identical(centroid[match(members$barcode,barcode)],members$centroid),
      identical(calls$central_barcode,centroids$central_barcode),
      all(calls$sum_counts==centroids$sum_counts))
    tm <- previous$timing
    tm$clustering_seconds <- unname(elapsed['elapsed'])
    tm$clustering_cpu_seconds <- unname(sum(elapsed[c('user.self','sys.self')]))
    mode <- 'Release rerun: every membership and centroid count matched development results'
  }
  result <- list(signature=signature,well=well,status='passed',mode=mode,
    input_hashes=hashes,timing=tm)
  saveRDS(result,cache)
  cat(well,'release verified:',tm$clustering_seconds,'elapsed seconds\n')
  result
}
checks <- parallel::mclapply(wells,check,mc.cores=workers,mc.preschedule=FALSE)
stopifnot(all(vapply(checks,function(x)is.list(x)&&identical(x$status,'passed'),logical(1))))
readr::write_csv(do.call(rbind,lapply(checks,`[[`,'timing')),file.path(kit,'results/release_clustering_times.csv'))
sources <- c('benchmark_release.R','load_release.R','build_release.R','provenance_helpers.R')
json_checks <- lapply(checks,function(x){x$input_hashes <- as.list(x$input_hashes);x})
jsonlite::write_json(list(status='passed',build=build,workers=workers,
  source_hashes=as.list(setNames(vapply(sources,function(f)sha(file.path(kit,f)),character(1)),sources)),
  all_memberships_and_centroid_counts_identical=TRUE,populations=json_checks,
  timing_sha256=sha(file.path(kit,'results/release_clustering_times.csv'))),
  file.path(kit,'results/release_validation.json'),pretty=TRUE,auto_unbox=TRUE)
cat('All six populations verified with the optimized release binary.\n')
