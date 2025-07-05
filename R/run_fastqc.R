#' Run FastQC on all fastq files in the sample table
#'
#' @param sample_table A data.frame or tibble with columns `R1` and `R2` pointing to FASTQ files
#' @param output_dir Directory where FastQC output will be written. Default is "fastQC"
#' 
#' @return Runs FastQC via system calls. Returns a character vector of commands run.
#' @export
run_fastqc <- function(sample_table, output_dir = "fastQC") {
  if (!all(c("R1", "R2") %in% colnames(sample_table))) {
    stop("sample_table must contain 'R1' and 'R2' columns.")
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  fastq_files <- unique(c(sample_table$R1, sample_table$R2))
  
  commands <- vapply(fastq_files, function(fq) {
    cmd <- sprintf("fastqc %s -o %s", fq, output_dir)
    message("Running: ", cmd)
    system(cmd)
    return(cmd)
  }, character(1))
  
  invisible(commands)
}
