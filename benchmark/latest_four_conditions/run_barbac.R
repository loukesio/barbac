# A fresh process for each method/distance/input; both exports are timed.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
.libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages(library(barbac))
cat(sprintf("BARBAC_BUILD_ID=%s\n", barbac:::barbac_build_id()))
elapsed <- system.time(result <- super_cluster2(
  args[[2]], method = args[[4]], distance = as.integer(args[[5]]),
  tie_break = "support", merge_ratio = 20, error_rate = 0.005,
  indel_model = args[[6]], use_design = FALSE, verbose = TRUE
))[["elapsed"]]
cat(sprintf("BARBAC_ALGO_SECONDS=%.6f\n", elapsed))
readr::write_csv(result[c("central_barcode", "sum_counts")],
                 file.path(args[[3]], "centroids.csv"))
readr::write_csv(data.frame(
  member = unlist(result$all_barcodes, use.names = FALSE),
  central_barcode = rep(result$central_barcode, lengths(result$all_barcodes)),
  member_count = unlist(result$all_counts, use.names = FALSE)
), file.path(args[[3]], "members.csv"))
