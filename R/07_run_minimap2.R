#' Map merged reads to reference using minimap2 and process with samtools
#'
#' @param merged_dir Directory where merged fastq files are stored (default: "merged")
#' @param reference Path to reference FASTA file
#' @param output_dir Directory for BAM files (default: "merged/bam")
#'
#' @return Vector of commands run
#' @export
run_minimap2 <- function(merged_dir = "merged", reference, output_dir = file.path(merged_dir, "bam")) {
  if (!file.exists(reference)) {
    stop("Reference file does not exist: ", reference)
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # \u2705 Match both _assembled.fastq and .assembled.fastq
  merged_files <- list.files(merged_dir, pattern = "([_.])assembled\\.fastq$", full.names = TRUE)
  if (length(merged_files) == 0) {
    stop("No merged files found with either _assembled.fastq or .assembled.fastq suffix")
  }
  
  commands <- vapply(merged_files, function(fq) {
    base <- tools::file_path_sans_ext(basename(fq))
    sam <- file.path(output_dir, paste0(base, ".sam"))
    bam <- file.path(output_dir, paste0(base, ".bam"))
    sorted_bam <- file.path(output_dir, paste0(base, "_sorted.bam"))
    
    cmds <- c(
      sprintf("minimap2 -a %s %s > %s",
              shQuote(reference), shQuote(fq), shQuote(sam)),
      sprintf("samtools view -S -b %s > %s", shQuote(sam), shQuote(bam)),
      sprintf("samtools sort %s -o %s", shQuote(bam), shQuote(sorted_bam)),
      sprintf("samtools index %s", shQuote(sorted_bam))
    )
    
    message("Running mapping pipeline for: ", fq)
    for (cmd in cmds) {
      status <- system(cmd)
      if (status != 0L) {
        stop("Mapping command failed with exit code ", status,
             ": ", cmd)
      }
    }
    
    return(paste(cmds, collapse = " && "))
  }, character(1))
  
  invisible(commands)
}
