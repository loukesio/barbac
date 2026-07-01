# Parameterised barbac runner.
# Usage: Rscript run_barbac.R <input.csv> <out_centroids.csv> [distance] [merge_ratio]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: run_barbac.R <input.csv> <out.csv> [distance] [merge_ratio]")
input_csv   <- args[[1]]
out_csv     <- args[[2]]
distance    <- if (length(args) >= 3) as.numeric(args[[3]]) else 3
merge_ratio <- if (length(args) >= 4) as.numeric(args[[4]]) else 20

suppressPackageStartupMessages({
  library(barbac)
  library(readr)
})

# Time just the super_cluster2() call, separately from R boot + package
# loading. The wrapper script measures the full wall clock; this timing
# isolates the algorithm itself so the small-dataset benchmark isn't
# dominated by ~4 s of R startup.
t0 <- Sys.time()
result <- super_cluster2(input_csv, distance = distance, merge_ratio = merge_ratio)
algo_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

out <- data.frame(
  central_barcode = result$central_barcode,
  sum_counts      = result$sum_counts
)
write_csv(out, out_csv)
cat(sprintf("Wrote %d centroids to %s\n", nrow(out), out_csv))
cat(sprintf("BARBAC_ALGO_SECONDS=%.4f\n", algo_s))
