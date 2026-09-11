#!/usr/bin/env Rscript
# How much of a clustering result rests on the arbitrary order of count ties?
#
# Barcodes that share a read count have no natural order, yet whichever is
# visited first is the one allowed to seed a cluster. super_cluster2() breaks
# ties deterministically so a result never depends on input row order, but the
# rule it uses is still one arbitrary choice among many. Re-running with
# tie_break = "hash" across several tie_seed values re-draws that choice and
# turns it into a measurable spread, which is the honest error bar to read the
# small between-method differences against.
#
# Writes one row per (condition, tie order) to tiebreak_seeds.csv.
#
# Usage:
#   Rscript run_tiebreak_sensitivity.R [results_dir] [n_seeds]

suppressPackageStartupMessages({
  library(barbac)
  library(readr)
  library(stringdist)
})

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else
  file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "indel_experiment", "results")
n_seeds     <- if (length(args) >= 2) as.integer(args[[2]]) else 15L

MAX_DIST    <- 3
MERGE_RATIO <- 20
CONDITIONS  <- c("sub_only", "low_indel",
                 "sub_only_structured", "low_indel_structured")

# Same definitions as the Python benchmark harness, so the numbers are directly
# comparable to the published tables:
#   FN  true barcodes absent from the centroid set
#   FP  centroids that are not true barcodes
#   WS  false positives lying within MAX_DIST of some true barcode
evaluate <- function(centroids, truth) {
  true_set <- unique(truth$barcode)
  cent_set <- unique(centroids$central_barcode)

  fn <- setdiff(true_set, cent_set)
  fp <- setdiff(cent_set, true_set)

  ws <- 0L
  if (length(fp) > 0) {
    # Chunked so the distance matrix stays small on the larger conditions.
    for (i in seq(1, length(fp), by = 200)) {
      chunk <- fp[i:min(i + 199, length(fp))]
      dm <- stringdist::stringdistmatrix(chunk, true_set, method = "lv",
                                         nthread = 1)
      ws <- ws + sum(apply(dm, 1, min) <= MAX_DIST)
    }
  }

  matched <- merge(centroids, truth,
                   by.x = "central_barcode", by.y = "barcode")
  r <- if (nrow(matched) >= 2)
    stats::cor(log10(matched$sum_counts), log10(matched$true_count)) else NA_real_

  list(n_centroids = nrow(centroids), n_true = length(true_set),
       fn = length(fn), fp = length(fp), ws = ws, pearson_r = r)
}

rows <- list()
for (cond in CONDITIONS) {
  input <- file.path(results_dir, cond, "input.csv")
  truth_path <- file.path(results_dir, cond, "true_counts.csv")
  if (!file.exists(input)) {
    message("skipping ", cond, " (no input.csv)")
    next
  }
  truth <- readr::read_csv(truth_path, show_col_types = FALSE)
  names(truth) <- c("barcode", "true_count")

  # tie_seed 0 with tie_break = "sequence" is the shipped default; the hash
  # seeds are equally valid alternative orders of the same tied barcodes.
  settings <- c(list(list(tie_break = "sequence", tie_seed = 0L)),
                lapply(seq_len(n_seeds),
                       function(s) list(tie_break = "hash", tie_seed = s)))

  for (cfg in settings) {
    t0 <- Sys.time()
    res <- super_cluster2(input, distance = MAX_DIST, merge_ratio = MERGE_RATIO,
                          tie_break = cfg$tie_break, tie_seed = cfg$tie_seed,
                          verbose = FALSE)
    algo_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    ev <- evaluate(data.frame(central_barcode = res$central_barcode,
                              sum_counts      = res$sum_counts), truth)
    rows[[length(rows) + 1]] <- data.frame(
      condition   = cond,
      tie_break   = cfg$tie_break,
      tie_seed    = cfg$tie_seed,
      is_default  = cfg$tie_break == "sequence",
      n_centroids = ev$n_centroids,
      n_true      = ev$n_true,
      fn = ev$fn, fp = ev$fp, ws = ev$ws,
      fn_pct = 100 * ev$fn / ev$n_true,
      fp_pct = 100 * ev$fp / ev$n_true,
      ws_pct = 100 * ev$ws / ev$n_true,
      pearson_r = ev$pearson_r,
      algo_s = algo_s,
      stringsAsFactors = FALSE
    )
    message(sprintf("%-22s %-9s seed=%-3d FN=%4d FP=%5d WS=%5d (%.2fs)",
                    cond, cfg$tie_break, cfg$tie_seed, ev$fn, ev$fp, ev$ws,
                    algo_s))
  }
}

out <- do.call(rbind, rows)
dest <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                   value = TRUE)[1])),
                  "tiebreak_seeds.csv")
readr::write_csv(out, dest)
message("\nwrote ", dest, "  (", nrow(out), " rows)")
