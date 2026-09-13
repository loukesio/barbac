#!/usr/bin/env Rscript
# One complete real-data FASTQ-to-plot run. Invoked by run.py after input checks.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
input_root <- normalizePath(args[1]); out <- normalizePath(args[2])
.libPaths(c(normalizePath(args[3]), .libPaths()))
clock <- proc.time()[['elapsed']]
stages <- list()
stage <- function(name, code) {
  start <- proc.time()[['elapsed']]
  cat(format(Sys.time(), tz='UTC', usetz=TRUE), name, '\n'); flush.console()
  value <- force(code)
  stages[[length(stages)+1L]] <<- data.frame(stage=name,
    elapsed_seconds=proc.time()[['elapsed']]-start)
  write.csv(do.call(rbind,stages),file.path(out,'stages.csv'),row.names=FALSE)
  value
}
stage('package_loading', {library(barbac); invisible(NULL)})
manifest <- read.delim(file.path(input_root,'samples.tsv'))
manifest <- manifest[manifest$well=='A3' & manifest$passage %in% c(2,4,6),]
stopifnot(nrow(manifest)==12, all(manifest$library_layout=='SINGLE'))
manifest$sample <- sprintf('A3_p%02d',manifest$passage)
reference <- file.path(input_root,'reference/cassette.fasta')
samples <- stage('input_staging', {
  dir.create(file.path(out,'inputs'))
  rows <- lapply(unique(manifest$sample), function(sample) {
    records <- manifest[manifest$sample==sample,]
    target <- file.path(out,'inputs',paste0(sample,'.fastq.gz'))
    dst <- file(target,'wb')
    tryCatch(for (accession in records$run_accession) {
      src <- file(file.path(input_root,'generated/fastq',paste0(accession,'.fastq.gz')),'rb')
      tryCatch(repeat {
        chunk <- readBin(src,'raw',n=4*1024^2)
        if (!length(chunk)) break
        writeBin(chunk,dst)
      }, finally=close(src))
    }, finally=close(dst))
    data.frame(sample=sample,R1=target,time=unique(records$passage)*6)
  })
  do.call(rbind,rows)
})
pipeline <- stage('fastq_to_bam_and_qc',
  run_cli_pipeline(samples,reference,file.path(out,'pipeline'),verbose=TRUE))
stopifnot(sum(pipeline$stats$mapped+pipeline$stats$unmapped)==sum(manifest$read_count),
          pipeline$multiqc_status=='completed', all(pipeline$command_timings$exit_status==0))
counts <- stage('barcode_extraction', {
  parts <- lapply(seq_len(nrow(samples)),function(i) {
    path <- barbac_xtr(pipeline$bam_files[[samples$sample[i]]],
      ref_name='Jasinska2020_barcode_cassette',start_pos=32,end_pos=40,
      flank_pattern='^([ACGT]{10,20})TATCTCGGTAG',min_count=1,verbose=FALSE,
      output_file=file.path(out,paste0(samples$sample[i],'-barcodes.csv')))
    x <- read.csv(path)
    stopifnot(sum(x$counts)==attr(path,'extraction_stats')$matched_alignments,nrow(x)>0)
    x$sample <- samples$sample[i]; x$time <- samples$time[i]; x
  })
  do.call(rbind,parts)
})
pooled <- stage('count_pooling', {
  x <- dplyr::summarise(dplyr::group_by(counts,barcode),counts=sum(counts),.groups='drop')
  stopifnot(sum(x$counts)==sum(counts$counts)); x
})
settings <- list(method='lv',distance=3,merge_ratio=20,error_rate=.005,
                 tie_break='support',indel_model='poisson')
clusters <- stage('native_clustering', do.call(super_cluster2,
  c(list(input_path=pooled,verbose=FALSE),settings)))
stats <- stage('cluster_statistics',cluster_stats(clusters,verbose=FALSE))
series <- stage('sample_assignment_and_export', {
  ids <- sprintf('L%06d',seq_len(nrow(clusters)))
  members <- data.frame(cluster_id=rep(ids,lengths(clusters$all_barcodes)),
    barcode=unlist(clusters$all_barcodes,use.names=FALSE))
  stopifnot(!anyDuplicated(members$barcode),setequal(members$barcode,pooled$barcode),
            sum(clusters$sum_counts)==sum(counts$counts))
  counts$cluster_id <- members$cluster_id[match(counts$barcode,members$barcode)]
  x <- dplyr::summarise(dplyr::group_by(counts,sample,time,cluster_id),counts=sum(counts),.groups='drop')
  stopifnot(identical(tapply(x$counts,x$sample,sum),tapply(counts$counts,counts$sample,sum)))
  readr::write_csv(x,file.path(out,'time_series.csv'))
  readr::write_csv(members,file.path(out,'memberships.csv'))
  readr::write_csv(stats,file.path(out,'cluster_stats.csv'))
  saveRDS(clusters,file.path(out,'clusters.rds'))
  x
})
stage('lineage_plot_and_export', {
  p <- barbac_ts_area(series,id_col='cluster_id',min_total_count=0,fill_missing='zero',
    palette='alger',x_breaks=samples$time,x_lab='Generation',
    y_lab='Fraction of extracted barcode reads',
    title='E. coli · low chloramphenicol, replicate 1')
  ggplot2::ggsave(file.path(out,'lineage_trajectories.pdf'),p,width=7,height=3.6,
                 device=grDevices::pdf,bg='white')
  ggplot2::ggsave(file.path(out,'lineage_trajectories.png'),p,width=7,height=3.6,dpi=300,bg='white')
  invisible(NULL)
})
result <- list(status='passed',package_version=as.character(packageVersion('barbac')),
  build_id=barbac:::barbac_build_id(),R=R.version.string,settings=settings,
  elapsed_seconds=proc.time()[['elapsed']]-clock,input_reads=sum(manifest$read_count),
  mapped_reads=sum(pipeline$stats$mapped),extracted_reads=sum(counts$counts),
  unique_sequences=nrow(pooled),clusters=nrow(clusters),samples=samples[c('sample','time')],
  sample_counts=as.list(tapply(counts$counts,counts$sample,sum)),
  no_count_loss=TRUE,all_lineages_plotted=TRUE,multiqc_status=pipeline$multiqc_status,
  scope='Raw local R1 reads; all four technical runs per sample; no quality filtering or UMI deduplication. Entire count matrix plotted, no abundance cutoff. Downloads, installation and input verification excluded; input concatenation included.')
jsonlite::write_json(result,file.path(out,'result.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
cat('COMPLETE:',result$input_reads,'reads;',result$elapsed_seconds,'seconds\n')
