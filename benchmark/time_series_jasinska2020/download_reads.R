#!/usr/bin/env Rscript
root <- normalizePath('benchmark/time_series_jasinska2020')
dir.create(file.path(root,'generated/fastq'),showWarnings=FALSE,recursive=TRUE)
x <- rbind(read.delim(file.path(root,'samples.tsv')),read.delim(file.path(root,'baseline_samples.tsv')))
args <- commandArgs(TRUE)
if(length(args)) x <- x[x$run_accession %in% args, ]
fetch <- function(i) {
  row <- x[i, ]; path <- file.path(root,'generated/fastq',paste0(row$run_accession,'.fastq.gz'))
  valid <- function() file.exists(path) && file.info(path)$size==row$fastq_bytes &&
    unname(tools::md5sum(path)) == row$fastq_md5
  if(!valid()) {
    status <- system2('curl',c('-fLsS','--retry','4','--connect-timeout','30',
      shQuote(paste0('https://',row$fastq_ftp)),'-o',shQuote(paste0(path,'.part'))))
    if(status!=0L) stop('Failed ',row$run_accession)
    if(file.info(paste0(path,'.part'))$size!=row$fastq_bytes ||
       unname(tools::md5sum(paste0(path,'.part'))) != row$fastq_md5) stop('Checksum failed ',row$run_accession)
    if(!file.rename(paste0(path,'.part'),path)) stop('Could not finalize ',path)
  }
  cat(row$run_accession,'verified\n')
  TRUE
}
out <- parallel::mclapply(seq_len(nrow(x)),fetch,mc.cores=4L,mc.preschedule=FALSE)
stopifnot(all(vapply(out,isTRUE,logical(1))))
