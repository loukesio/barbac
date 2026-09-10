#!/usr/bin/env Rscript
# Before/after real-data check for exact search pruning; never tune assignments.
devtools::load_all(quiet=TRUE)
root <- normalizePath('benchmark/time_series_jasinska2020')
out <- file.path(root,'generated/search_check');dir.create(out,recursive=TRUE,showWarnings=FALSE)
phase <- commandArgs(TRUE)[1]
if(phase=='baseline') {
  x <- readr::read_csv(file.path(root,'generated/clustering/A3/input.csv'),show_col_types=FALSE)
  set.seed(20200910)
  ids <- c(seq_len(10000),sample(10001:nrow(x),2000))
  input <- x[ids, ]
  saveRDS(input,file.path(out,'input.rds'))
} else input <- readRDS(file.path(out,'input.rds'))
canonical <- function(x) {
  z <- data.frame(barcode=unlist(x$all_barcodes,use.names=FALSE),
    centroid=rep(x$central_barcode,lengths(x$all_barcodes)),
    count=unlist(x$all_counts,use.names=FALSE))
  z <- z[order(z$barcode), ];rownames(z)<-NULL;z
}
tm <- system.time(x <- barbac::super_cluster2(input,method='lv',distance=3,tie_break='support',
  merge_ratio=20,error_rate=.005,indel_model='poisson',use_design=FALSE,verbose=FALSE))
record <- list(build=barbac:::barbac_build_id(),timing=unclass(tm),assignments=canonical(x),
  clusters=x[c('central_barcode','sum_counts')])
saveRDS(record,file.path(out,paste0(phase,'.rds')))
if(phase!='baseline') {
  baseline <- readRDS(file.path(out,'baseline.rds'))
  stopifnot(identical(record$assignments,baseline$assignments),
    identical(record$clusters$central_barcode,baseline$clusters$central_barcode),
    identical(record$clusters$sum_counts,baseline$clusters$sum_counts))
  jsonlite::write_json(list(status='passed',sequences=nrow(input),
    memberships_identical=TRUE,centroids_and_counts_identical=TRUE,
    before_build=baseline$build,after_build=record$build,
    before_cpu_seconds=sum(baseline$timing[1:2]),after_cpu_seconds=sum(record$timing[1:2]),
    before_elapsed_seconds=baseline$timing[3],after_elapsed_seconds=record$timing[3]),
    file.path(root,'sources/search_validation.json'),pretty=TRUE,auto_unbox=TRUE)
}
print(tm)
