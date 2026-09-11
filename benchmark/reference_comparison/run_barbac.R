# Process time includes startup and both exports; core time includes CSV reading.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
.libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages(library(barbac))
cat(sprintf('BARBAC_BUILD_ID=%s\n', barbac:::barbac_build_id()))
elapsed <- system.time(result <- super_cluster2(
  args[[2]], method = args[[4]], tie_break = args[[5]],
  distance = 3, merge_ratio = 20, error_rate = 0.005, verbose = TRUE
))[["elapsed"]]
cat(sprintf('BARBAC_ALGO_SECONDS=%.6f\n', elapsed))
readr::write_csv(result[c('central_barcode', 'sum_counts')],
                 file.path(args[[3]], 'centroids.csv'))
members <- data.frame(
  member = unlist(result$all_barcodes, use.names = FALSE),
  central_barcode = rep(result$central_barcode, lengths(result$all_barcodes)),
  member_count = unlist(result$all_counts, use.names = FALSE)
)
readr::write_csv(members, file.path(args[[3]], 'members.csv'))
