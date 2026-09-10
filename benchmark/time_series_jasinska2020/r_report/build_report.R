#!/usr/bin/env Rscript
repo <- normalizePath('.')
kit <- file.path(repo,'benchmark/time_series_jasinska2020');here <- file.path(kit,'r_report')
source(file.path(kit,'load_release.R'));invisible(load_release(kit))
source(file.path(kit,'provenance_helpers.R'))
source(file.path(here,'report_helpers.R'))
validation <- jsonlite::read_json(file.path(kit,'results/validation.json'),simplifyVector=TRUE)
stopifnot(validation$status=='passed')
assert_hash_map(validation$source_hashes,kit)
assert_hash_map(validation$result_hashes,file.path(kit,'results'))
report <- make_report_assets(read_report(kit),kit,here)
saveRDS(report,file.path(kit,'generated/report_data.rds'))
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(here,'session_info.txt'))
setwd(here)
status <- system2('quarto',c('render','report.qmd','--to','html','--output','report.html','--self-contained','--quiet'))
stopifnot(status==0L)
files <- c('report.qmd','report_helpers.R','build_report.R','report.html','session_info.txt')
jsonlite::write_json(list(status='numerical_checks_passed',browser_validation='pending',
  hashes=as.list(setNames(vapply(files,function(f)digest::digest(file=f,algo='sha256'),character(1)),files)),
  no_competitor_runs=TRUE,publication_comparison='Sample summaries and published dominant-barcode final frequencies',
  publication_normalization_resolved=FALSE,
  explicit_comparison_denominators=c('All input reads','Extracted barcode reads'),
  samples=75,biological_populations=6,published_dominant_barcodes=nrow(report$top)),
  'validation.json',pretty=TRUE,auto_unbox=TRUE)
