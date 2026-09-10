blue <- '#2C4A63'
gold <- '#B07D22'
report_gt <- function(data,title,subtitle=NULL) {
  gt::gt(data) |>
    gt::tab_header(title=title,subtitle=subtitle) |>
    gt::tab_options(table.font.names='Arial',table.font.size=13,
      heading.title.font.size=18,table.width=gt::pct(100),data_row.padding=gt::px(7))
}
report_widget <- function(plot,height=480) {
  w <- plotly::ggplotly(plot,tooltip=c('text','x','y'),height=height,width=NULL) |>
    plotly::layout(autosize=TRUE) |>
    plotly::config(displaylogo=FALSE,responsive=TRUE)
  w$width <- '100%';w$sizingPolicy$knitr$figure <- FALSE;w$sizingPolicy$defaultWidth <- '100%'
  w
}
report_explorer <- function(data) {
  DT::datatable(data,rownames=FALSE,extensions='Buttons',filter='top',
    colnames=c('Population','Barcode','Passage','Barbac / input reads (%)',
      'Barbac / extracted reads (%)','Published final (%)'),
    options=list(dom='Bfrtip',pageLength=10,scrollX=TRUE,deferRender=TRUE,
      buttons=list(list(extend='csv',text='Download filtered CSV',
        exportOptions=list(orthogonal='export'))))) |>
    DT::formatRound(columns=4:6,digits=4)
}
read_report <- function(kit,include_runtime=TRUE) {
  source(file.path(kit,'provenance_helpers.R'),local=TRUE)
  results <- file.path(kit,'results')
  provenance <- jsonlite::read_json(file.path(results,'provenance.json'),simplifyVector=TRUE)
  stopifnot(provenance$status=='complete')
  assert_hash_map(provenance$source_hashes,kit)
  assert_hash_map(provenance$package_sources,file.path(kit,'../..'))
  r <- list(qc=readr::read_csv(file.path(results,'extraction_summary.csv'),show_col_types=FALSE),
    diversity=readr::read_csv(file.path(results,'diversity.csv'),show_col_types=FALSE),
    comparison=readr::read_csv(file.path(results,'publication_sample_comparison.csv'),show_col_types=FALSE),
    top=readr::read_csv(file.path(results,'publication_top_barcodes.csv'),show_col_types=FALSE),
    provenance=provenance)
  if(include_runtime) {
    r$times <- readr::read_csv(file.path(results,'release_clustering_times.csv'),show_col_types=FALSE)
    release <- jsonlite::read_json(file.path(results,'release_validation.json'),simplifyVector=TRUE)
    stopifnot(release$status=='passed',release$all_memberships_and_centroid_counts_identical,
      release$timing_sha256==digest::digest(file=file.path(results,'release_clustering_times.csv'),algo='sha256'))
    assert_hash_map(release$source_hashes,kit)
    assert_hash_map(release$build$source_hashes,file.path(kit,'../..'))
    r$release <- release
  }
  r$qc$treatment <- sub('Chloramphenicol 1 microgram/mL','Chloramphenicol',r$qc$treatment,fixed=TRUE)
  r$diversity$treatment <- sub('Chloramphenicol 1 microgram/mL','Chloramphenicol',r$diversity$treatment,fixed=TRUE)
  paths <- list.files(file.path(kit,'generated/clustering'),pattern='time_series.rds$',recursive=TRUE,full.names=TRUE)
  r$series <- lapply(paths,readRDS);names(r$series) <- vapply(r$series,`[[`,character(1),'population')
  stopifnot(length(r$series)==6,nrow(r$comparison)==75,
    all(r$comparison$input_reads_barbac==r$comparison$input_reads_publication),
    all(r$qc$quality_passed==r$qc$mapped+r$qc$unmapped),
    all(r$qc$extracted_reads<=r$qc$mapped))
  source(file.path(kit,'frequency_audit.R'),local=TRUE)
  frequency <- compare_frequency_denominators(r$top,r$comparison)
  r$top <- frequency$top;r$final_summary <- frequency$summary;r$denominator_audit <- frequency$audit
  r
}
make_report_assets <- function(report,kit,here) {
  dir.create(file.path(here,'figures'),showWarnings=FALSE)
  cache_dir <- file.path(kit,'generated/report_diagnostics')
  dir.create(cache_dir,showWarnings=FALSE)
  # Diagnostics for the first and last measured timepoint of every population.
  q <- report$qc[report$qc$passage>0, ]
  chosen <- do.call(rbind,lapply(split(q,q$population),function(x)x[x$passage %in% range(x$passage), ]))
  diagnostics <- lapply(seq_len(nrow(chosen)),function(i) {
    row <- chosen[i, ]
    data_file <- file.path(kit,'generated/samples',row$sample,'barcodes.csv')
    figures <- file.path(here,'figures',paste0(row$sample,'_',1:3,'.png'))
    sha <- function(f)digest::digest(file=f,algo='sha256')
    signature <- digest::digest(c(sha(data_file),sha(file.path(here,'report_helpers.R')),
      sha(file.path(kit,'../../R/10_barbac_xtr.R'))),algo='sha256')
    cache <- file.path(cache_dir,paste0(row$sample,'.rds'))
    if(file.exists(cache) && all(file.exists(figures))) {
      old <- readRDS(cache)
      if(identical(signature,old$signature) && identical(vapply(figures,sha,character(1)),old$figure_hashes))
        return(old$diagnostic)
    }
    data <- readr::read_csv(data_file,show_col_types=FALSE)
    details <- barbac::barbac_xtr.stats(data,c(14,16),fill_color=blue,verbose=FALSE,
      panel_labels=TRUE,return_details=TRUE)
    for(j in seq_along(details$plots)) {
      p <- details$plots[[j]] + barbac::theme_barbac(base_size=12,family='Arial') +
        ggplot2::labs(tag=LETTERS[j])
      if(j==3) p <- p+ggplot2::labs(x='Base-composition entropy (bits)',title='Entropy within each barcode')
      ggplot2::ggsave(file.path(here,'figures',paste0(row$sample,'_',j,'.png')),p,
        width=6,height=4,dpi=130,device=ragg::agg_png)
    }
    diagnostic <- list(sample=row$sample,population=row$population,passage=row$passage,summary=details$length_summary)
    saveRDS(list(signature=signature,figure_hashes=vapply(figures,sha,character(1)),diagnostic=diagnostic),cache)
    diagnostic
  })
  report$diagnostics <- diagnostics
  source(file.path(here,'composition_helpers.R'),local=TRUE)
  # One shared mapping ensures a barcode keeps its colour across all six panels.
  ids <- sort(unique(unlist(lapply(report$series,function(s)rownames(s$count_matrix)))))
  colours <- setNames(viridisLite::plasma(length(ids)),ids)
  readr::write_csv(data.frame(barcode=ids,colour=unname(colours)),
    file.path(here,'barcode_colours.csv.gz'))
  # Separate R processes are safe for macOS font rendering; forked workers can
  # abort in Objective-C font initialization and return NULL without an R error.
  workers <- parallel::makePSOCKcluster(3L)
  on.exit(parallel::stopCluster(workers),add=TRUE)
  invisible(parallel::clusterCall(workers,function(kit,here) {
    source(file.path(kit,'load_release.R'));invisible(load_release(kit))
    source(file.path(here,'composition_helpers.R'),local=globalenv())
    NULL
  },kit,here))
  report$composition <- parallel::parLapply(workers,report$series,
    function(s,colours,here)make_composition(s,colours,here),colours,here)
  stopifnot(length(report$composition)==length(report$series),
    all(vapply(report$composition,function(x)is.list(x) &&
      isTRUE(x$every_frequency_and_band_width_verified) &&
      identical(x$lineages,x$polygon_count),logical(1))))
  explorer <- list()
  for(pop in names(report$series)) {
    s <- report$series[[pop]];m <- s$count_matrix
    top <- s$top
    for(i in seq_len(nrow(top))) explorer[[length(explorer)+1L]] <- data.frame(
      Population=pop,Barcode=top$barcode[i],Passage=s$passages,
      Barbac_frequency_input_percent=100*s$top_frequencies[i,],
      Barbac_frequency_extracted_percent=100*s$top_frequencies[i,]*s$totals/colSums(m),
      Published_final_percent=ifelse(s$passages==max(s$passages),100*top$published_final[i],NA_real_))
  }
  report$explorer <- do.call(rbind,explorer)
  jsonlite::write_json(list(status='passed',populations=report$composition,
    shared_barcode_colours=TRUE,colour_map_sha256=digest::digest(
      file=file.path(here,'barcode_colours.csv.gz'),algo='sha256')),
    file.path(here,'composition_validation.json'),pretty=TRUE,auto_unbox=TRUE)
  # Keep large matrices in separate reproducible count files, not the HTML cache.
  report$series <- NULL
  report
}
