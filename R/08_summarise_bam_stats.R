#' Summarise mapped and unmapped read counts from sorted BAM files
#'
#' @param bam_dir Directory containing sorted BAM files (default: "merged/bam")
#'
#' @return A data.frame with sample names, mapped and unmapped counts
#' @export
summarise_bam_stats <- function(bam_dir = "merged/bam") {
  bam_files <- list.files(bam_dir, pattern = "_sorted\\.bam$", full.names = TRUE)
  if (length(bam_files) == 0) stop("No _sorted.bam files found in ", bam_dir)

  stats <- lapply(bam_files, function(bam) {
    sample <- gsub("_sorted\\.bam$", "", basename(bam))
    mapped_out <- system2("samtools", c("view", "-c", "-F", "4", shQuote(bam)),
                          stdout = TRUE, stderr = TRUE)
    unmapped_out <- system2("samtools", c("view", "-c", "-f", "4", shQuote(bam)),
                            stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(mapped_out, "status")) || !is.null(attr(unmapped_out, "status"))) {
      stop("samtools failed while summarising ", bam, ".")
    }
    mapped <- as.integer(mapped_out[[1L]])
    unmapped <- as.integer(unmapped_out[[1L]])
    tibble::tibble(sample = sample, mapped = mapped, unmapped = unmapped)
  })

  result <- do.call(rbind, stats)
  return(result)
}

#' Per-sample bar plot of mapping counts
#'
#' Companion visualisation for [summarise_bam_stats()]. Draws a bar plot of
#' mapped vs unmapped read counts across samples so outlier samples are
#' visible at a glance.
#'
#' @param stats Either a data.frame with columns `sample`, `mapped`,
#'   `unmapped` (as returned by [summarise_bam_stats()]) or a path to a
#'   directory of sorted BAM files (in which case `summarise_bam_stats()`
#'   is called first).
#' @param position Character. Bar layout: `"stack"` (default) or `"dodge"`.
#' @param mapped_color Character. Fill colour for mapped reads.
#'   Default `"#2C4A63"` (barbac primary).
#' @param unmapped_color Character. Fill colour for unmapped reads.
#'   Default `"#B0B7BF"`.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip labs scale_fill_manual
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
plot_bam_stats <- function(stats,
                           position       = c("stack", "dodge"),
                           mapped_color   = "#2C4A63",
                           unmapped_color = "#B0B7BF") {
  position <- match.arg(position)

  if (is.character(stats) && length(stats) == 1L && dir.exists(stats)) {
    stats <- summarise_bam_stats(stats)
  }
  if (!is.data.frame(stats)) {
    stop("`stats` must be a data.frame or a directory containing sorted BAM files.")
  }
  missing_cols <- setdiff(c("sample", "mapped", "unmapped"), names(stats))
  if (length(missing_cols)) {
    stop("`stats` is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  long <- tidyr::pivot_longer(stats,
                              cols      = c("mapped", "unmapped"),
                              names_to  = "status",
                              values_to = "reads")
  long$status <- factor(long$status, levels = c("mapped", "unmapped"))

  ggplot2::ggplot(long, ggplot2::aes(x    = .data$sample,
                                     y    = .data$reads,
                                     fill = .data$status)) +
    ggplot2::geom_bar(stat = "identity", position = position) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(mapped   = mapped_color,
                                          unmapped = unmapped_color)) +
    ggplot2::labs(x     = NULL,
                  y     = "Reads",
                  fill  = NULL,
                  title = "Per-sample mapping counts")
}
