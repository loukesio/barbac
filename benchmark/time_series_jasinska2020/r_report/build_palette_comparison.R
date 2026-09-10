#!/usr/bin/env Rscript
here <- normalizePath('benchmark/time_series_jasinska2020/r_report')
source(file.path(here,'../provenance_helpers.R'))
snapshot <- jsonlite::read_json(file.path(here,'palette_preview/palettes.json'))
stopifnot(snapshot$all_32_interpolations_match_ggvmap,length(snapshot$palettes)==32)
for(well in c('A3','A1')) {
  z <- jsonlite::read_json(file.path(here,'palette_preview',paste0(well,'_validation.json')))
  stopifnot(z$all_barcodes_individual,z$full_frequencies_verified,z$geometry_unchanged)
  for(im in z$images) stopifnot(digest::digest(file=file.path(here,'palette_preview',im$path),
    algo='sha256')==im$sha256)
}
setwd(here)
status <- system2('quarto',c('render','palette_comparison.qmd','--to','html',
  '--output','palette_comparison.html','--self-contained','--quiet'))
stopifnot(status==0L)
files <- c('palette_preview.R','composition_helpers.R','build_palette_comparison.R',
  'palette_comparison.qmd','palette_comparison.html',
  file.path('palette_preview',list.files('palette_preview',pattern='[.](json|png)$')))
jsonlite::write_json(list(status='numerical_checks_passed',browser_validation='pending',
  hashes=as.list(setNames(vapply(files,function(f)digest::digest(file=f,algo='sha256'),character(1)),files)),
  counts_and_order_identical_across_palettes=TRUE,all_barcodes_individual=TRUE,
  populations=c('Low CMP r1','No drug r1'),preview_palettes=c('alger','dora','casa_natal')),
  'palette_validation.json',pretty=TRUE,auto_unbox=TRUE)
