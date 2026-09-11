#!/usr/bin/env Rscript
kit <- normalizePath('benchmark/time_series_jasinska2020')
manifest <- rbind(read.delim(file.path(kit,'samples.tsv')),
                  read.delim(file.path(kit,'baseline_samples.tsv')))
stopifnot(nrow(manifest)==312,!anyDuplicated(manifest$run_accession))
files <- list.files(file.path(kit,'generated/samples'),pattern='receipt[.]json$',recursive=TRUE,full.names=TRUE)
stopifnot(length(files)==78)
receipts <- lapply(files,jsonlite::read_json,simplifyVector=TRUE)
sha <- function(f)digest::digest(file=f,algo='sha256')
seen <- character()
for(r in receipts) {
  rows <- manifest[match(r$runs,manifest$run_accession),]
  stopifnot(nrow(rows)==4,!anyNA(rows$run_accession),
    sum(rows$read_count)==r$input_reads,
    r$reference_sha256==sha(file.path(kit,'reference/cassette.fasta')),
    r$script_sha256==sha(file.path(kit,'process_samples.R')),
    r$barcodes_sha256==sha(file.path(kit,'generated/samples',r$sample,'barcodes.csv')),
    identical(unname(r$input_md5),unname(rows$fastq_md5)),
    r$mapped+r$unmapped==r$quality_passed,r$quality_passed<=r$input_reads,
    r$extracted_reads<=r$mapped)
  for(i in seq_len(nrow(rows))) {
    run <- rows$run_accession[i]
    zip <- file.path(kit,'generated/samples',r$sample,'fastqc',paste0(run,'_fastqc.zip'))
    con <- unz(zip,paste0(run,'_fastqc/fastqc_data.txt'))
    lines <- readLines(con,warn=FALSE);close(con)
    n <- as.numeric(sub('^Total Sequences\t','',lines[startsWith(lines,'Total Sequences\t')]))
    stopifnot(length(n)==1,n==rows$read_count[i])
  }
  seen <- c(seen,r$runs)
}
stopifnot(!anyDuplicated(seen),setequal(seen,manifest$run_accession))
out <- file.path(kit,'results');dir.create(out,showWarnings=FALSE)
jsonlite::write_json(receipts,file.path(out,'processing_receipts.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
jsonlite::write_json(list(status='passed',samples=78,FASTQs=312,
  input_reads=sum(manifest$read_count),fastqc_totals_reconciled=312,
  sample_receipt_hashes=as.list(setNames(vapply(files,sha,character(1)),
    vapply(receipts,`[[`,character(1),'sample')))),file.path(out,'processing_validation.json'),pretty=TRUE,auto_unbox=TRUE)
cat('All 312 independent FastQC totals and 78 processing receipts reconciled.\n')
