#' Run MultiQC to summarize FastQC reports
#'
#' @param input_dir The directory where FastQC output is located
#' @param output_dir The directory where MultiQC summary should be saved
#'
#' @return Runs MultiQC and returns the command used
#' @export
run_multiqc <- function(input_dir = "fastQC", output_dir = "multiQC") {
  if (.barbac_which("multiqc") == "") {
    stop("MultiQC not found in system PATH.")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cmd <- sprintf("multiqc %s -o %s", shQuote(input_dir), shQuote(output_dir))
  message("Running: ", cmd)
  status <- .barbac_system(cmd)
  if (status != 0L) {
    stop("MultiQC failed with exit code ", status, ".")
  }
  
  invisible(cmd)
}
