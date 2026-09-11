#!/usr/bin/env Rscript
# Public source inventory and deterministic cohort selection, before analysis.
root <- normalizePath('benchmark/time_series_jasinska2020')
fetch <- function(url, path) {
  if (!file.exists(path)) {
    status <- system2('curl', c('-fLsS', '--retry', '3', shQuote(url), '-o', shQuote(path)))
    if (status != 0L) stop('Download failed: ', url)
  }
  invisible(path)
}
library(xml2)
library(readxl)
aliases <- read.delim(file.path(root, 'sources/constant_aliases.tsv'))
runs <- read.delim(file.path(root, 'sources/PRJNA592529_runs.tsv'))
x <- merge(aliases, runs, by='run_accession')
x$well <- sub('.*well-([^_]+).*', '\\1', x$sample_alias)
x$passage <- as.integer(sub('.*passage-([0-9]+).*', '\\1', x$sample_alias))
x$subsample <- as.integer(sub('.*subsample-([0-9]+).*', '\\1', x$sample_alias))
conditions <- lapply(sort(unique(x$well)), function(well) {
  row <- x[x$well == well, ][1, ]
  path <- file.path(root, 'sources', paste0(row$sample_accession,'.xml'))
  fetch(paste0('https://www.ebi.ac.uk/ena/browser/api/xml/',row$sample_accession),path)
  doc <- read_xml(path)
  value <- xml_text(xml_find_first(doc,"//SAMPLE_ATTRIBUTE[TAG='Drug condition and replicate']/VALUE"))
  data.frame(well=well, population=value, example_sample=row$sample_accession)
})
conditions <- do.call(rbind,conditions)
write.csv(conditions,file.path(root,'sources/well_conditions.csv'),row.names=FALSE)
x <- merge(x,conditions,by='well')
selected <- x[grepl('^(Low CMP|No drug) r[123]$',x$population), ]
selected <- selected[order(selected$population,selected$passage,selected$subsample), ]
selected$treatment <- ifelse(grepl('^Low CMP',selected$population),'Chloramphenicol 1 microgram/mL','No antibiotic')
selected$replicate <- as.integer(sub('.*r([123])$','\\1',selected$population))
write.table(selected,file.path(root,'samples.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
baseline <- x[x$well == 'all', ]
baseline$population <- 'Shared initial population'
baseline$treatment <- 'Shared baseline'
baseline$replicate <- as.integer(sub('.*sample-([0-9]+)_subsample.*','\\1',baseline$sample_alias))
write.table(baseline,file.path(root,'baseline_samples.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
book <- file.path(root,'sources/supplementary_tables.xlsx')
for (s in c('Supplementary Table 1c','Supplementary Table 4b')) {
  tab <- suppressMessages(read_excel(book,sheet=s,col_names=FALSE))
  write.csv(tab,file.path(root,'sources',paste0(gsub(' ','_',s),'.csv')),row.names=FALSE)
}
print(conditions,row.names=FALSE)
cat('\nSelected:',nrow(selected),'runs;',sum(selected$read_count),'reads;',sum(selected$fastq_bytes)/1e9,'GB compressed\n')
print(table(selected$population,selected$passage))
