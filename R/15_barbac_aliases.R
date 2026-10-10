# Preferred public names for the core workflow functions.
#
# Each alias forwards every argument unchanged to the original function. The
# original names remain exported and fully supported: existing scripts,
# benchmark receipts, Studio and tests continue to call them. The validated
# clustering sources (R/11_super_cluster2.R, src/clustering.cpp and the Rcpp
# export files) are not modified by these aliases.

#' Cluster barcodes (preferred name for `super_cluster2()`)
#'
#' `barbac_cluster()` is the public name for the native abundance-aware
#' centroid clustering engine. It forwards all arguments to
#' [super_cluster2()] and returns the same object. Both names are supported.
#'
#' @param ... Arguments passed to [super_cluster2()], such as `input_path`,
#'   `distance`, `method`, `merge_ratio`, `error_rate`, `tie_break` and
#'   `indel_model`.
#' @return The tibble returned by [super_cluster2()]: one row per cluster with
#'   `cluster_id`, `central_barcode`, `all_barcodes`, `all_counts` and
#'   `sum_counts`.
#' @seealso [super_cluster2()] for the full argument documentation,
#'   [cluster_stats()] for post-clustering diagnostics.
#' @examples
#' input <- data.frame(
#'   barcode = c("AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAT"),
#'   counts = c(1000, 10)
#' )
#' barbac_cluster(input, distance = 3, verbose = FALSE)
#' @export
barbac_cluster <- function(...) {
  super_cluster2(...)
}

#' Extract barcodes from BAM files (preferred name for `barbac_xtr()`)
#'
#' `barbac_extract()` forwards all arguments to [barbac_xtr()] and returns
#' the same result. Both names are supported.
#'
#' @param ... Arguments passed to [barbac_xtr()], such as `bam_file`,
#'   reference coordinates or a flank pattern, and `output_file`.
#' @return The value returned by [barbac_xtr()].
#' @seealso [barbac_xtr()] for the full argument documentation,
#'   [barbac_extract_stats()] for summaries of the extracted barcodes.
#' @export
barbac_extract <- function(...) {
  barbac_xtr(...)
}

#' Summarise extracted barcodes (preferred name for `barbac_xtr.stats()`)
#'
#' `barbac_extract_stats()` forwards all arguments to [barbac_xtr.stats()]
#' and returns the same result. Both names are supported.
#'
#' @param ... Arguments passed to [barbac_xtr.stats()], such as `file` and
#'   `barcode_length`.
#' @return The value returned by [barbac_xtr.stats()].
#' @seealso [barbac_xtr.stats()] for the full argument documentation.
#' @export
barbac_extract_stats <- function(...) {
  barbac_xtr.stats(...)
}

#' Map reads to a reference (preferred name for `run_cli_pipeline()`)
#'
#' `barbac_map()` forwards all arguments to [run_cli_pipeline()], which runs
#' FastQC, optional PEAR merging, minimap2 and samtools and writes sorted,
#' indexed BAM files with mapping statistics. Both names are supported.
#'
#' @param ... Arguments passed to [run_cli_pipeline()], such as the sample
#'   table, reference and output directory.
#' @return The invisible list returned by [run_cli_pipeline()].
#' @seealso [run_cli_pipeline()] for the full argument documentation.
#' @export
barbac_map <- function(...) {
  run_cli_pipeline(...)
}
