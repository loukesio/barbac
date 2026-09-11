# One isolated R process per timed run, including CSV loading in algorithm time.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop('Expected library, input, output, method, tie_break')
.libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages(library(barbac))
elapsed <- system.time(result <- super_cluster2(
  args[[2]], method = args[[4]], tie_break = args[[5]], verbose = FALSE
))[["elapsed"]]
write.csv(result[c('central_barcode', 'sum_counts')], args[[3]], row.names = FALSE)
cat(sprintf('BARBAC_ALGO_SECONDS=%.6f\n', elapsed))
cat(sprintf('BARBAC_BUILD_ID=%s\n', barbac:::barbac_build_id()))
