#' Run FastQC on all FASTQ files listed in a samples.csv file
#'
#' This function reads a `samples.csv` file from a given path, extracts `R1` and optionally `R2` FASTQ paths,
#' and runs FastQC on all unique files. The output is saved under `results/fastQC/` in the same directory as the CSV.
#'
#' @param sample_csv Path to the `samples.csv` file.
#'
#' @return A character vector of FastQC system commands executed.
#' @importFrom utils read.csv
#' @export
run_fastqc <- function(sample_csv) {
  # Validate the CSV path
  if (!file.exists(sample_csv)) {
    stop("\u274C Cannot find the samples.csv file at: ", sample_csv)
  }
  
  # Extract base directory where samples.csv is located
  base_dir <- dirname(sample_csv)
  
  # Read sample table
  sample_table <- read.csv(sample_csv)
  
  # Check for required column
  if (!("R1" %in% colnames(sample_table))) {
    stop("\u274C samples.csv must contain at least an 'R1' column.")
  }
  
  # Combine R1 and R2 if available
  fastq_files <- sample_table$R1
  if ("R2" %in% colnames(sample_table)) {
    fastq_files <- unique(c(fastq_files, sample_table$R2))
  } else {
    message("\u2139 Only 'R1' column found \u2014 single-end mode.")
  }
  
  # Check if FastQC is available
  if (Sys.which("fastqc") == "") {
    warning("\u26A0 FastQC not found in system PATH. Please install it or patch your PATH.")
    return(invisible(character(0)))
  }
  
  # Create output directory inside base_dir/results/fastQC
  output_dir <- file.path(base_dir, "results", "fastQC")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Run FastQC
  commands <- vapply(fastq_files, function(fq) {
    cmd <- sprintf("fastqc %s -o %s", fq, output_dir)
    message("\u25B6 Running: ", cmd)
    system(cmd)
    cmd
  }, character(1))
  
  invisible(commands)
}

#barbac::run_fastqc("path/to/samples.csv")
