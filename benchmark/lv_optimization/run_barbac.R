# Process time includes startup and both exports; core time includes CSV reading.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 5:7)
options_extra <- if (length(args) >= 6L) list(indel_model = args[[6]]) else list()
rate <- if (length(args) >= 7L) as.numeric(args[[7]]) else 0.005
.libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages(library(barbac))
cat(sprintf('BARBAC_BUILD_ID=%s\n', barbac:::barbac_build_id()))
elapsed <- system.time(result <- do.call(super_cluster2, c(list(
  args[[2]], method = args[[4]], tie_break = args[[5]],
  distance = 3, merge_ratio = 20, error_rate = rate, verbose = TRUE
), options_extra)))[["elapsed"]]
cat(sprintf('BARBAC_ALGO_SECONDS=%.6f\n', elapsed))
readr::write_csv(result[c('central_barcode', 'sum_counts')],
                 file.path(args[[3]], 'centroids.csv'))
members <- data.frame(
  member = unlist(result$all_barcodes, use.names = FALSE),
  central_barcode = rep(result$central_barcode, lengths(result$all_barcodes)),
  member_count = unlist(result$all_counts, use.names = FALSE)
)
readr::write_csv(members, file.path(args[[3]], 'members.csv'))
