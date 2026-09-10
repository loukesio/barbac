#!/usr/bin/env Rscript
# Independent reconciliation of exported counts, receipts and author tables.
kit <- normalizePath('benchmark/time_series_jasinska2020')
source(file.path(kit,'analysis_helpers.R'))
source(file.path(kit,'provenance_helpers.R'))
source(file.path(kit,'r_report/report_helpers.R'))
r <- read_report(kit,include_runtime=FALSE)
manifest <- read.delim(file.path(kit,'samples.tsv'))
baseline <- read.delim(file.path(kit,'baseline_samples.tsv'))
stopifnot(nrow(manifest)==300,nrow(baseline)==12,
  sum(manifest$read_count)==sum(r$qc$input_reads[r$qc$passage>0]),
  sum(baseline$read_count)==sum(r$qc$input_reads[r$qc$passage==0]),
  nrow(r$top)==120,all(table(r$top$population)==20))
author <- read_author_tables(kit)
expected <- author$top[author$top$population %in% r$top$population,]
key <- function(x) paste(x$population,x$barcode)
stopifnot(setequal(key(expected),key(r$top)))
idx <- match(key(r$top),key(expected))
stopifnot(identical(r$top$published_final,expected$published_final[idx]))
raw_baseline <- lapply(1:3,function(i) readr::read_csv(file.path(kit,'generated/samples',
  paste0('initial_',i),'barcodes.csv'),show_col_types=FALSE))
checks <- list()
for(pop in names(r$series)) {
  s <- r$series[[pop]];well <- s$manifest$well[1]
  members <- readr::read_csv(file.path(kit,'generated/clustering',well,'members.csv'),show_col_types=FALSE)
  exported <- readr::read_csv(file.path(kit,'results',paste0(well,'_barcode_counts.csv.gz')),show_col_types=FALSE)
  stopifnot(identical(exported$barcode,rownames(s$count_matrix)),
    identical(dim(as.matrix(exported[-1])),dim(s$count_matrix)),
    all(unname(as.matrix(exported[-1]))==unname(s$count_matrix)),
    !anyDuplicated(members$barcode))
  inputs <- c(list(do.call(rbind,raw_baseline)),lapply(s$manifest$sample,function(id)
    readr::read_csv(file.path(kit,'generated/samples',id,'barcodes.csv'),show_col_types=FALSE)))
  for(j in seq_along(inputs)) {
    z <- inputs[[j]];centroid <- members$centroid[match(z$barcode,members$barcode)]
    stopifnot(!anyNA(centroid))
    expected_counts <- rowsum(z$counts,centroid,reorder=FALSE)
    observed <- s$count_matrix[match(rownames(expected_counts),rownames(s$count_matrix)),j]
    stopifnot(all(observed==as.numeric(expected_counts)),
      sum(s$count_matrix[,j])==sum(z$counts))
    d <- r$diversity[r$diversity$population==pop & r$diversity$passage==s$passages[j],]
    positive <- s$count_matrix[,j];positive <- positive[positive>0]
    p <- positive/sum(positive)
    stopifnot(nrow(d)==1,d$richness==length(positive),
      isTRUE(all.equal(d$assigned_shannon_effective,exp(-sum(p*log(p))))),
      isTRUE(all.equal(d$assigned_inverse_dominance,1/max(p))))
  }
  top <- r$top[r$top$population==pop,]
  matched <- match(top$barcode,rownames(s$count_matrix))
  freq <- rep(0,nrow(top));freq[!is.na(matched)] <-
    s$count_matrix[matched[!is.na(matched)],ncol(s$count_matrix)]/tail(s$totals,1)
  stopifnot(isTRUE(all.equal(freq,top$barbac_final)),
    isTRUE(all.equal(freq*tail(s$totals,1)/sum(s$count_matrix[,ncol(s$count_matrix)]),top$barbac_final_extracted)))
  checks[[pop]] <- list(population=pop,timepoints=ncol(s$count_matrix),
    clusters=nrow(s$count_matrix),all_count_cells_reconciled=TRUE,
    shared_baseline_extracted_reads=sum(s$count_matrix[,1]))
}
readr::write_csv(r$top,file.path(kit,'results/publication_frequency_comparison.csv'))
readr::write_csv(r$final_summary,file.path(kit,'results/publication_frequency_summary.csv'))
readr::write_csv(r$denominator_audit,file.path(kit,'results/publication_denominator_audit.csv'))
cat('All count cells, diversity values and publication frequencies reconciled; verifying release timing provenance.\n')
deadline <- Sys.time()+7200
while(!file.exists(file.path(kit,'results/release_validation.json'))) {
  if(Sys.time()>deadline)stop('Release verification is not complete')
  Sys.sleep(10)
}
release_report <- read_report(kit)
release <- jsonlite::read_json(file.path(kit,'results/release_validation.json'))
stopifnot(length(release$populations)==6)
for(check in release$populations) {
  stopifnot(check$status=='passed')
  assert_hash_map(check$input_hashes,file.path(kit,'generated/clustering',check$well))
}
files <- list.files(file.path(kit,'results'),full.names=TRUE,pattern='\\.(csv|gz|json)$')
files <- files[basename(files)!='validation.json']
jsonlite::write_json(list(status='passed',sample_read_totals_match_publication=75,
  published_dominant_barcode_entries=120,unique_input_reads=sum(r$qc$input_reads),
  populations=checks,full_count_matrices_reconciled=TRUE,
  author_final_frequency_normalization='unresolved; both explicit denominators reported',
  source_hashes=as.list(setNames(vapply(c('validate_results.R','frequency_audit.R','r_report/report_helpers.R','provenance_helpers.R'),
    function(f)digest::digest(file=file.path(kit,f),algo='sha256'),character(1)),
    c('validate_results.R','frequency_audit.R','r_report/report_helpers.R','provenance_helpers.R'))),
  result_hashes=as.list(setNames(vapply(files,function(f)digest::digest(file=f,algo='sha256'),character(1)),basename(files)))),
  file.path(kit,'results/validation.json'),pretty=TRUE,auto_unbox=TRUE)
cat('All sample totals, count-matrix cells, diversity values and 120 published barcode frequencies reconciled.\n')
