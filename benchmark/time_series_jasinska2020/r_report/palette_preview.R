#!/usr/bin/env Rscript
# Optional full-data palette previews; does not alter the main report palette.
kit <- normalizePath('benchmark/time_series_jasinska2020')
here <- file.path(kit,'r_report');out <- file.path(here,'palette_preview')
dir.create(out,showWarnings=FALSE)
# Use the current plotting API. The historical installed release and its
# clustering timing receipts remain frozen at their recorded source hashes.
repo <- normalizePath(file.path(kit,'../..'))
pkgload::load_all(repo,quiet=TRUE,export_all=FALSE)
source(file.path(here,'composition_helpers.R'))
source(file.path(kit,'provenance_helpers.R'))
palettes <- barbac::barbac_palettes()
args <- commandArgs(TRUE);well <- if(length(args))args[1] else 'swatches'
chosen <- c('alger','dora','casa_natal')
expand_palette <- barbac::barbac_palette
if(well=='swatches') {
  stopifnot(length(palettes)==32,
    identical(expand_palette('Casa Natal',100),expand_palette('casa_natal',100)))
  # Compare against the archived colour data used for the existing previews.
  # This requires no installation of the package from which they originated.
  original <- jsonlite::read_json(file.path(out,'original_palettes.json'),simplifyVector=TRUE)
  stopifnot(identical(palettes,original$palettes))
  for(name in names(palettes)) stopifnot(identical(expand_palette(name,100),
    grDevices::colorRampPalette(original$palettes[[name]])(100)))
  shared <- readr::read_csv(file.path(here,'barcode_colours.csv.gz'),show_col_types=FALSE)$barcode
  for(name in chosen) stopifnot(identical(expand_palette(name,length(shared)),
    grDevices::colorRampPalette(original$palettes[[name]])(length(shared))))
  d <- do.call(rbind,lapply(names(palettes),function(name)data.frame(
    palette=name,x=seq_len(200),colour=expand_palette(name,200))))
  d$palette <- factor(d$palette,levels=rev(names(palettes)))
  p <- ggplot2::ggplot(d,ggplot2::aes(x,palette,fill=colour))+
    ggplot2::geom_raster()+ggplot2::scale_fill_identity()+
    ggplot2::scale_x_continuous(expand=c(0,0),breaks=NULL)+
    ggplot2::labs(title='All 32 LTC palettes built into barbac',
      subtitle='Each row interpolates the palette end to end; these are colour options, not measured values.',
      x=NULL,y=NULL)+ggplot2::theme_minimal(base_size=12,base_family='Arial')+
    ggplot2::theme(panel.grid=ggplot2::element_blank(),plot.margin=ggplot2::margin(12,25,12,14))
  ggplot2::ggsave(file.path(out,'all_palettes.png'),p,width=10,height=12,dpi=160,device=ragg::agg_png)
  files <- c('R/12_barbac_ts_area.R','R/13_theme_barbac.R','R/14_barbac_palettes.R',
    'NAMESPACE','DESCRIPTION')
  jsonlite::write_json(list(source='barbac::barbac_palettes()',
    package_version=as.character(utils::packageVersion('barbac')),
    source_hashes=as.list(setNames(vapply(files,function(f)
      digest::digest(file=file.path(repo,f),algo='sha256'),character(1)),files)),
    palettes=palettes,all_32_interpolations_match_original=TRUE,
    shared_preview_colour_vectors_identical=TRUE,shared_barcodes=length(shared),
    preview_geometry_and_images_retained=TRUE,
    curated_fills=TRUE),file.path(out,'palettes.json'),pretty=TRUE,auto_unbox=TRUE)
  cat('All 32 palette interpolations verified and swatches saved.\n')
  quit(status=0)
}
stopifnot(well %in% c('A3','A1'))
s <- readRDS(file.path(kit,'generated/clustering',well,'time_series.rds'))
m <- s$count_matrix
shared <- readr::read_csv(file.path(here,'barcode_colours.csv.gz'),show_col_types=FALSE)$barcode
ids <- sort(rownames(m));stopifnot(all(ids %in% shared))
input <- data.frame(barcode=rep(rownames(m),ncol(m)),
  time=rep(s$passages,each=nrow(m)),counts=as.vector(m))
message('Preparing full-data palette comparison: ',s$population)
p <- barbac::barbac_ts_area(input,min_total_count=0,fill_missing='zero',
  palette=chosen[1],show_legend=FALSE,
  x_breaks=c(0,6,12,18,24,30),x_lab='Passage',y_lab='Barcode frequency among extracted reads',
  title=if(well=='A3')'Chloramphenicol · replicate 1' else 'No antibiotic · replicate 1',
  theme=barbac::theme_barbac(base_size=12,family='Arial')+
    ggplot2::theme(plot.margin=ggplot2::margin(10,25,12,14)))
stack <- composition_stack(p)
expected <- sweep(m[stack$ids,,drop=FALSE],2,colSums(m),'/')
stopifnot(isTRUE(all.equal(unname(stack$freq),unname(expected),tolerance=1e-12)))
receipt <- list(population=s$population,well=well,lineages=nrow(m),
  all_barcodes_individual=TRUE,min_total_count=0,
  denominator='All extracted barcode reads',full_frequencies_verified=TRUE,
  geometry_sha256=digest::digest(stack,algo='sha256'),images=list())
for(name in chosen) {
  colours <- setNames(expand_palette(name,length(shared)),shared)
  q <- suppressMessages(p+ggplot2::scale_fill_manual(values=colours[ids]))+
    ggplot2::labs(subtitle=paste(format(nrow(m),big.mark=',',trim=TRUE),'barcodes ·',name))
  stopifnot(identical(unname(q$scales$get_scales('fill')$palette(length(ids))),unname(colours[ids])))
  relative <- paste0(well,'_',name,'.png')
  ggplot2::ggsave(file.path(out,relative),draw_composition_stack(q,stack),
    width=9,height=6,dpi=180,device=ragg::agg_png)
  receipt$images[[name]] <- list(path=relative,
    sha256=digest::digest(file=file.path(out,relative),algo='sha256'))
  message('Saved ',relative)
}
receipt$geometry_unchanged <- identical(receipt$geometry_sha256,digest::digest(stack,algo='sha256'))
stopifnot(receipt$geometry_unchanged)
jsonlite::write_json(receipt,file.path(out,paste0(well,'_validation.json')),pretty=TRUE,auto_unbox=TRUE)
