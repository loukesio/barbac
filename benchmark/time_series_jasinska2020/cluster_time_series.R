#!/usr/bin/env Rscript
root <- normalizePath('benchmark/time_series_jasinska2020')
source(file.path(root,'load_release.R'));release_build <- load_release(root)
source(file.path(root,'analysis_helpers.R'))
sha <- function(path) digest::digest(file=path,algo='sha256')
out <- file.path(root,'results');dir.create(out,showWarnings=FALSE)
manifest <- read.delim(file.path(root,'samples.tsv'))
requested <- commandArgs(TRUE)
if(length(requested)) manifest <- manifest[manifest$well %in% requested, ]
stopifnot(nrow(manifest)>0)
manifest$sample <- sprintf('%s_p%02d',manifest$well,manifest$passage)
sample_ids <- unique(c(manifest$sample,paste0('initial_',1:3)))
receipts <- lapply(sample_ids,function(id) {
  file <- file.path(root,'generated/samples',id,'receipt.json')
  deadline <- Sys.time()+7200
  while(!file.exists(file)) {if(Sys.time()>deadline)stop('Missing ',id);Sys.sleep(10)}
  r <- jsonlite::read_json(file,simplifyVector=TRUE)
  stopifnot(r$barcodes_sha256==sha(file.path(root,'generated/samples',id,'barcodes.csv')))
  r
}); names(receipts) <- sample_ids
raw <- lapply(sample_ids,function(id) readr::read_csv(file.path(root,'generated/samples',id,'barcodes.csv'),show_col_types=FALSE))
names(raw) <- sample_ids
baseline <- collapse_counts(raw[paste0('initial_',1:3)])
baseline_total <- sum(vapply(receipts[paste0('initial_',1:3)],`[[`,numeric(1),'input_reads'))
author <- read_author_tables(root)
saveRDS(author,file.path(root,'generated/author_tables.rds'))
summaries <- comparisons <- times <- list()
for(pop in unique(manifest$population)) {
  m <- unique(manifest[manifest$population==pop,c('sample','well','passage','population','treatment','replicate')])
  m <- m[order(m$passage), ]
  inputs <- c(list(baseline=baseline),raw[m$sample])
  pooled <- collapse_counts(inputs)
  cdir <- file.path(root,'generated/clustering',m$well[1]);dir.create(cdir,recursive=TRUE,showWarnings=FALSE)
  readr::write_csv(pooled,file.path(cdir,'input.csv'))
  cache <- file.path(cdir,'completed.rds')
  signature <- digest::digest(list(input=sha(file.path(cdir,'input.csv')),
    source=sha(file.path(root,'../../src/clustering.cpp')),
    wrapper=sha(file.path(root,'../../R/11_super_cluster2.R')),
    settings=c('lv','3','support','20','.005','poisson','no-design')),algo='sha256')
  if(file.exists(cache)) {
    previous <- readRDS(cache)
    stopifnot(previous$signature==signature)
    summaries[[length(summaries)+1L]] <- previous$summary
    comparisons[[pop]] <- previous$comparison
    times[[pop]] <- previous$timing
    next
  }
  local_summary <- list()
  timing <- system.time(calls <- barbac::super_cluster2(pooled,barcode_col='barcode',counts_col='counts',
    method='lv',distance=3,tie_break='support',merge_ratio=20,error_rate=.005,
    indel_model='poisson',use_design=FALSE,verbose=TRUE))
  seconds <- timing[['elapsed']]
  cpu_seconds <- sum(timing[c('user.self','sys.self')])
  members <- data.frame(barcode=unlist(calls$all_barcodes,use.names=FALSE),
    centroid=rep(calls$central_barcode,lengths(calls$all_barcodes)))
  stopifnot(!anyDuplicated(members$barcode),setequal(members$barcode,pooled$barcode),
            sum(calls$sum_counts)==sum(pooled$counts))
  readr::write_csv(members,file.path(cdir,'members.csv'))
  readr::write_csv(calls[c('central_barcode','sum_counts')],file.path(cdir,'centroids.csv'))
  ids <- calls$central_barcode
  count_matrix <- matrix(0,nrow=length(ids),ncol=length(inputs),dimnames=list(ids,names(inputs)))
  totals <- c(baseline_total,vapply(receipts[m$sample],`[[`,numeric(1),'input_reads'))
  passages <- c(0,m$passage)
  for(j in seq_along(inputs)) {
    input <- inputs[[j]]
    centroid <- members$centroid[match(input$barcode,members$barcode)]
    stopifnot(!anyNA(centroid))
    cnt <- tapply(input$counts,centroid,sum)
    count_matrix[match(names(cnt),ids),j] <- cnt
    stopifnot(sum(count_matrix[,j])==sum(input$counts))
    d <- diversity(count_matrix[,j],totals[j])
    local_summary[[length(local_summary)+1L]] <- data.frame(population=pop,treatment=m$treatment[1],
      replicate=m$replicate[1],passage=passages[j],generation=6*passages[j],
      input_reads=totals[j],assigned_reads=sum(cnt),as.list(d),
      as.list(barbac::cluster_stats(data.frame(sum_counts=unname(cnt)),verbose=FALSE)))
  }
  colnames(count_matrix) <- paste0('passage_',passages)
  readr::write_csv(data.frame(barcode=ids,count_matrix,check.names=FALSE),
    file.path(out,paste0(m$well[1],'_barcode_counts.csv.gz')))
  top <- author$top[author$top$population==pop, ]
  idx <- match(top$barcode,ids)
  freqs <- matrix(0,nrow=nrow(top),ncol=ncol(count_matrix))
  freqs[!is.na(idx), ] <- sweep(count_matrix[idx[!is.na(idx)],,drop=FALSE],2,totals,'/')
  top$barbac_final <- freqs[,ncol(freqs)]
  top$barbac_mean_including_baseline <- rowMeans(freqs)
  top$barbac_mean_sampled_only <- rowMeans(freqs[,-1,drop=FALSE])
  comparisons[[pop]] <- top
  times[[pop]] <- data.frame(population=pop,unique_input_sequences=nrow(pooled),clusters=nrow(calls),
    clustering_seconds=seconds,clustering_cpu_seconds=cpu_seconds,
    input_read_observations=sum(pooled$counts))
  saveRDS(list(population=pop,manifest=m,passages=passages,totals=totals,
    count_matrix=count_matrix,top=top,top_frequencies=freqs),file.path(cdir,'time_series.rds'))
  local_summary <- do.call(rbind,local_summary)
  summaries[[length(summaries)+1L]] <- local_summary
  saveRDS(list(signature=signature,summary=local_summary,comparison=top,timing=times[[pop]],
    binary_sha256=release_build$binary_sha256,build_mode='optimized release'),cache)
  cat(pop,'complete:',nrow(calls),'clusters;',seconds,'seconds\n')
}
if(length(requested)) {
  cat('Selected populations cached; full aggregation is performed by running without arguments.\n')
  quit(status=0L)
}
readr::write_csv(do.call(rbind,summaries),file.path(out,'diversity.csv'))
readr::write_csv(do.call(rbind,comparisons),file.path(out,'publication_top_barcodes.csv'))
readr::write_csv(do.call(rbind,times),file.path(out,'clustering_times.csv'))
fields <- c('sample','population','treatment','replicate','passage','generation','input_reads',
  'quality_passed','mapped','unmapped','extracted_reads','raw_sequences','quality_seconds',
  'mapping_seconds','extraction_seconds','fastqc_seconds','total_seconds')
qc <- do.call(rbind,lapply(receipts,function(r)as.data.frame(r[fields])))
readr::write_csv(qc,file.path(out,'extraction_summary.csv'))
comparison <- merge(qc,author$samples,by=c('population','passage'),suffixes=c('_barbac','_publication'))
stopifnot(nrow(comparison)==length(unique(manifest$sample)),
          all(comparison$input_reads_barbac==comparison$input_reads_publication))
readr::write_csv(comparison,file.path(out,'publication_sample_comparison.csv'))
jsonlite::write_json(list(status='complete',method='barbac LV + Poisson',build=barbac:::barbac_build_id(),distance=3,
  tie_break='support',merge_ratio=20,error_rate=.005,use_design=FALSE,
  baseline='Three shared initial samples pooled, reused once within each population clustering',
  primary_frequency_denominator='All FASTQ reads in the sample for read accounting; author final-frequency normalization is unresolved. The report also compares frequencies among extracted barcode reads.',
  reference='Published 288-base masked cassette; positions 11–25 are the barcode',
  package_sources=as.list(setNames(vapply(c('src/clustering.cpp','R/11_super_cluster2.R','R/10_barbac_xtr.R','R/10a_extract_flanked_bam.R'),
    function(f)sha(file.path(root,'../..',f)),character(1)),
    c('src/clustering.cpp','R/11_super_cluster2.R','R/10_barbac_xtr.R','R/10a_extract_flanked_bam.R'))),
  source_hashes=as.list(setNames(vapply(c('cluster_time_series.R','process_samples.R','analysis_helpers.R',
    'samples.tsv','baseline_samples.tsv','reference/cassette.fasta','sources/supplementary_tables.xlsx'),
    function(f)sha(file.path(root,f)),character(1)),c('cluster_time_series.R','process_samples.R','analysis_helpers.R',
    'samples.tsv','baseline_samples.tsv','reference/cassette.fasta','sources/supplementary_tables.xlsx')))),
  file.path(out,'provenance.json'),pretty=TRUE,auto_unbox=TRUE)
