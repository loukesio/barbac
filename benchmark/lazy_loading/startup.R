args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[[1]], .libPaths()))
elapsed <- system.time(suppressPackageStartupMessages(library(barbac)))[['elapsed']]
cat(sprintf('LOAD_SECONDS=%.6f\n', elapsed))
cat('GENOMIC_ALIGNMENTS_LOADED=', 'GenomicAlignments' %in% loadedNamespaces(), '\n', sep = '')
cat('BUILD_ID=', barbac:::barbac_build_id(), '\n', sep = '')
