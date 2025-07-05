#' Run MultiQC to summarize FastQC reports
#'
#' @param input_dir The directory where FastQC output is located
#' @param output_dir The directory where MultiQC summary should be saved
#'
#' @return Runs MultiQC and returns the command used
#' @export
run_multiqc <- function(input_dir = "fastQC", output_dir = "multiQC") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cmd <- sprintf("multiqc %s -o %s", input_dir, output_dir)
  message("Running: ", cmd)
  system(cmd)
  
  invisible(cmd)
}
