#' Post-clustering QC summary
#'
#' Reports high-level diagnostics from [super_cluster2()] output to help
#' sanity-check `distance` and `merge_ratio` choices before running
#' downstream time-series analysis. A very high singleton fraction usually
#' means `distance` was too tight; a dominant single cluster with a small
#' `top1_frac` << expected can indicate over-aggressive merging.
#'
#' @param clusters Tibble returned by [super_cluster2()], or a path to a
#'   CSV of that shape. Must contain a `sum_counts` column.
#' @param top_n Integer vector. Rank thresholds for cumulative-abundance
#'   share. Default `c(1, 10, 100)`.
#' @param singleton_max Integer. Maximum `sum_counts` treated as a
#'   singleton. Default `1L`.
#' @param verbose Logical. Print the summary to console. Default `TRUE`.
#'
#' @return A one-row tibble with columns:
#'   \itemize{
#'     \item `n_clusters`
#'     \item `total_reads`
#'     \item `n_singletons`, `singleton_frac`
#'     \item `largest_cluster_size`, `largest_cluster_frac`
#'     \item `topK_frac` for each `K` in `top_n`
#'   }
#'
#' @examples
#' \dontrun{
#' result <- super_cluster2("sample1_barcodes.csv", distance = 3)
#' cluster_stats(result)
#' }
#'
#' @export
#' @importFrom tibble as_tibble
#' @importFrom readr read_csv
cluster_stats <- function(clusters,
                          top_n         = c(1, 10, 100),
                          singleton_max = 1L,
                          verbose       = TRUE) {
  if (is.character(clusters) && length(clusters) == 1L) {
    if (!file.exists(clusters)) stop("File not found: ", clusters)
    clusters <- readr::read_csv(clusters, show_col_types = FALSE)
  }
  if (!is.data.frame(clusters)) {
    stop("`clusters` must be a data.frame/tibble or a CSV path.")
  }
  if (!"sum_counts" %in% names(clusters)) {
    stop("Input must contain a `sum_counts` column (super_cluster2 output).")
  }

  counts      <- sort(as.numeric(clusters$sum_counts), decreasing = TRUE)
  n_clusters  <- length(counts)
  total_reads <- sum(counts)
  if (n_clusters == 0L || total_reads == 0) {
    stop("No clusters / zero reads to summarise.")
  }

  n_singletons   <- sum(counts <= singleton_max)
  singleton_frac <- n_singletons / n_clusters
  largest_size   <- counts[1]
  largest_frac   <- largest_size / total_reads

  top_fracs <- vapply(top_n, function(k) {
    k <- min(k, n_clusters)
    sum(counts[seq_len(k)]) / total_reads
  }, numeric(1))
  names(top_fracs) <- paste0("top", top_n, "_frac")

  out <- tibble::as_tibble(c(
    list(
      n_clusters           = n_clusters,
      total_reads          = total_reads,
      n_singletons         = n_singletons,
      singleton_frac       = singleton_frac,
      largest_cluster_size = largest_size,
      largest_cluster_frac = largest_frac
    ),
    as.list(top_fracs)
  ))

  if (verbose) {
    fmt_pct <- function(x) sprintf("%.2f%%", 100 * x)
    message("Clusters:            ", format(n_clusters,   big.mark = ","))
    message("Total reads:         ", format(total_reads,  big.mark = ","))
    message("Singletons (<= ", singleton_max, "): ",
            format(n_singletons, big.mark = ","),
            " (", fmt_pct(singleton_frac), ")")
    message("Largest cluster:     ", format(largest_size, big.mark = ","),
            " (", fmt_pct(largest_frac), " of reads)")
    for (i in seq_along(top_n)) {
      message(sprintf("Top-%d abundance:     %s", top_n[i], fmt_pct(top_fracs[i])))
    }
  }

  out
}
