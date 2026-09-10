#!/usr/bin/env Rscript
# Single-end FASTQ -> quality filter -> reference mapping -> barbac extraction.
devtools::load_all(quiet=TRUE)
barbac::use_barbac_env()
root <- normalizePath('benchmark/time_series_jasinska2020')
args <- commandArgs(TRUE)
x <- rbind(read.delim(file.path(root,'samples.tsv')),read.delim(file.path(root,'baseline_samples.tsv')))
x$sample <- ifelse(x$well=='all',paste0('initial_',x$replicate),sprintf('%s_p%02d',x$well,x$passage))
groups <- split(x,x$sample)
groups <- groups[unique(x$sample)]
if(length(args)) groups <- groups[args]
stopifnot(length(groups)>0,!any(vapply(groups,is.null,logical(1))))
ref <- file.path(root,'reference/cassette.fasta')
sha <- function(path) digest::digest(file=path,algo='sha256')
process_sample <- function(rows) {
  sample <- rows$sample[1]
  out <- file.path(root,'generated/samples',sample)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  receipt <- file.path(out,'receipt.json')
  if(file.exists(receipt)) {
    old <- jsonlite::read_json(receipt,simplifyVector=TRUE)
    stopifnot(old$reference_sha256==sha(ref),old$script_sha256==sha(file.path(root,'process_samples.R')),
              old$barcodes_sha256==sha(file.path(out,'barcodes.csv')))
    return(old)
  }
  paths <- file.path(root,'generated/fastq',paste0(rows$run_accession,'.fastq.gz'))
  deadline <- Sys.time()+3600
  while(!all(file.exists(paths))) {
    if(Sys.time()>deadline) stop('Required downloads are incomplete for ',sample)
    Sys.sleep(5)
  }
  stopifnot(all(file.exists(paths)),identical(unname(tools::md5sum(paths)),rows$fastq_md5))
  begin <- proc.time()[['elapsed']]
  tempfq <- file.path(out,'quality_passed.fastq')
  if(file.exists(tempfq)) unlink(tempfq)
  total <- passed <- 0
  for(i in seq_along(paths)) {
    dna <- Biostrings::readDNAStringSet(paths[i],format='fastq',with.qualities=TRUE)
    stopifnot(length(dna)==rows$read_count[i])
    bad <- Biostrings::letterFrequency(S4Vectors::mcols(dna)$qualities,
      paste0(intToUtf8(33:42,multiple=TRUE),collapse=''),OR='|')[,1]
    keep <- bad==0L
    total <- total+length(dna); passed <- passed+sum(keep)
    Biostrings::writeXStringSet(dna[keep],tempfq,format='fastq',append=TRUE)
    rm(dna,bad,keep); gc(FALSE)
  }
  quality_seconds <- proc.time()[['elapsed']]-begin
  bam <- file.path(out,'mapped_sorted.bam')
  cmd <- paste('minimap2 -a -x sr -k 9 -w 5 -m 10 -s 10 -n 1 --secondary=no -t 2',
    shQuote(ref),shQuote(tempfq),'2>',shQuote(file.path(out,'minimap2.log')),
    '| samtools view -u -F 2304 - | samtools sort -@ 1 -m 256M -o',shQuote(bam),'-')
  tm <- proc.time()[['elapsed']]
  status <- system2('bash',c('-o','pipefail','-c',shQuote(cmd)))
  stopifnot(status==0L,system2('samtools',c('quickcheck',shQuote(bam)))==0L,
            system2('samtools',c('index',shQuote(bam)))==0L)
  mapping_seconds <- proc.time()[['elapsed']]-tm
  stats <- barbac::summarise_bam_stats(out)
  stopifnot(stats$mapped+stats$unmapped==passed)
  tm <- proc.time()[['elapsed']]
  extracted <- barbac::barbac_xtr(bam,ref_name='Jasinska2020_barcode_cassette',
    start_pos=32,end_pos=40,flank_pattern='^([ACGT]{10,20})TATCTCGGTAG',
    output_file=file.path(out,'barcodes.csv'),min_count=1,verbose=FALSE)
  data <- readr::read_csv(extracted,show_col_types=FALSE)
  extraction_seconds <- proc.time()[['elapsed']]-tm
  stopifnot(sum(data$counts)==attr(extracted,'extraction_stats')$matched_alignments,
            !anyDuplicated(data$barcode),all(data$barcode_length>=10 & data$barcode_length<=20))
  qdir <- file.path(out,'fastqc'); dir.create(qdir,showWarnings=FALSE)
  tm <- proc.time()[['elapsed']]
  status <- system2('fastqc',c('--quiet','--threads','2','--outdir',shQuote(qdir),shQuote(paths)),
    stdout=file.path(out,'fastqc.log'),stderr=file.path(out,'fastqc.err'))
  stopifnot(status==0L)
  fastqc_seconds <- proc.time()[['elapsed']]-tm
  result <- list(sample=sample,population=rows$population[1],well=rows$well[1],
    treatment=rows$treatment[1],replicate=rows$replicate[1],passage=rows$passage[1],
    generation=6*rows$passage[1],runs=rows$run_accession,input_reads=total,
    quality_passed=passed,mapped=stats$mapped,unmapped=stats$unmapped,
    extracted_reads=sum(data$counts),raw_sequences=nrow(data),
    quality_seconds=quality_seconds,mapping_seconds=mapping_seconds,
    extraction_seconds=extraction_seconds,fastqc_seconds=fastqc_seconds,
    total_seconds=proc.time()[['elapsed']]-begin,extraction=attr(extracted,'extraction_stats'),
    reference_sha256=sha(ref),script_sha256=sha(file.path(root,'process_samples.R')),
    barcodes_sha256=sha(extracted),input_md5=setNames(rows$fastq_md5,rows$run_accession))
  jsonlite::write_json(result,receipt,pretty=TRUE,auto_unbox=TRUE,digits=15)
  unlink(tempfq)
  cat(sample,'complete:',total,'reads;',sum(data$counts),'extracted\n')
  result
}
results <- parallel::mclapply(groups,process_sample,mc.cores=2L,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
