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
    mapped <- as.integer(system2("samtools", c("view", "-c", "-F", "4", bam), stdout = TRUE))
    unmapped <- as.integer(system2("samtools", c("view", "-c", "-f", "4", bam), stdout = TRUE))
    tibble::tibble(sample = sample, mapped = mapped, unmapped = unmapped)
  })
  
  result <- do.call(rbind, stats)
  return(result)
}
