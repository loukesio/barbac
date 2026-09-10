#!/usr/bin/env Rscript
# May run while clustering: these diagnostics depend only on completed extraction.
devtools::load_all(quiet=TRUE)
kit <- normalizePath('benchmark/time_series_jasinska2020');here <- file.path(kit,'r_report')
source(file.path(here,'report_helpers.R'));source(file.path(kit,'analysis_helpers.R'))
receipts <- jsonlite::read_json(file.path(kit,'results/processing_receipts.json'),simplifyVector=FALSE)
qc <- do.call(rbind,lapply(receipts,function(x)as.data.frame(x[c('sample','population','passage')],stringsAsFactors=FALSE)))
top <- read_author_tables(kit)$top
invisible(make_report_assets(list(qc=qc,top=top[top$population %in% qc$population,],series=list()),kit,here))
cat('Extraction diagnostic figures cached.\n')
