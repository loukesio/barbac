#' barbac: End-to-end DNA barcode lineage analysis
#'
#' `barbac` provides an R-native pipeline for tracking DNA barcode
#' lineages from raw sequencing reads to time-resolved lineage
#' trajectories. It bundles command-line pre-processing (FastQC, PEAR,
#' minimap2, samtools) via [run_cli_pipeline()] with a native R
#' clustering routine [super_cluster2()] and time-series visualisation
#' [barbac_ts_area()].
#'
#' @section Core clustering:
#'   [super_cluster2()] implements abundance-ranked greedy Levenshtein
#'   clustering with a distance-aware count-ratio merge guard, backed by
#'   a bit-parallel 64-bit edit-distance kernel written in C++. On the
#'   Johnson et al. (2023) benchmark it reaches statistical parity with
#'   Shepherd at roughly 1.7-fold the speed; on datasets with indel
#'   errors it substantially outperforms Shepherd on false-positive and
#'   wrong-sequence rates (see the package vignette).
#'
#' @section Pipeline:
#'   The full FASTQ-to-lineage workflow is orchestrated by
#'   [run_cli_pipeline()], which chains [run_fastqc()], [run_pear_merge()],
#'   [run_minimap2()], [summarise_bam_stats()], and optionally
#'   [run_multiqc()]. External tools are provisioned into a conda
#'   environment by [configure_environment()] and activated for the
#'   current R session by [use_barbac_env()].
#'
#' @section Visualisation:
#'   [barbac_ts_area()] renders stacked-area plots of lineage
#'   trajectories over time. It accepts long-format or Bartender-wide
#'   input, supports late-appearing lineages via grid completion, and
#'   fills missing (barcode, timepoint) cells with either zero or a
#'   small epsilon.
#'
#' @references
#'   Theodosiou L, Farr AD, Rainey PB (2026). *barbac: A versatile
#'   tool for analysing DNA barcode sequences.* Bioinformatics (in
#'   preparation).
#'
#'   Johnson MS, Venkataram S, Kryazhimskiy S (2023). *Best Practices
#'   in Designing, Sequencing, and Identifying Random DNA Barcodes.*
#'   J Mol Evol 91:263-280.
#'
#'   Tavakolian N, Cochran JR, Levy SF (2022). *Shepherd: accurate
#'   clustering for correcting DNA barcode errors.* Bioinformatics
#'   38:3710-3716.
#'
#' @seealso
#'   [super_cluster2()], [barbac_ts_area()], [barbac_xtr()],
#'   [run_cli_pipeline()]
#'
#' @examples
#' \dontrun{
#' # Cluster a barcode count table
#' result <- super_cluster2("reads.csv", distance = 3, merge_ratio = 20)
#'
#' # Visualise a time series of clustered barcode frequencies
#' barbac_ts_area(time_series_df, min_total_count = 10)
#' }
#'
#' @docType package
#' @name barbac-package
#' @aliases barbac
#' @useDynLib barbac, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

# Silence R CMD check "no visible binding for global variable" notes
# that come from tidyverse non-standard evaluation (column names used
# unquoted inside dplyr/tidyr/ggplot2 verbs).
utils::globalVariables(c(
  "barcode", "barcodes", "bin", "counts", "entropy", "reads", "x"
))
