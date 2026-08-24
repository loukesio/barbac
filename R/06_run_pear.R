#' Merge paired-end reads using PEAR
#'
#' This function reads a `samples.csv` file from the given path and merges reads using PEAR.
#' If only `R1` is available, the sample is skipped.
#' Output is saved to `results/merged/` inside the same folder.
#'
#' @param sample_csv Full path to `samples.csv`
#'
#' @return Character vector of PEAR commands run.
#' @importFrom utils read.csv
#' @export
run_pear_merge <- function(sample_csv) {
  # Validate path
  if (!file.exists(sample_csv)) {
    stop("\u274C samples.csv not found at: ", sample_csv)
  }
  
  # Read sample table and extract base directory
  sample_table <- read.csv(sample_csv)
  base_dir <- dirname(sample_csv)
  
  # Check required columns
  if (!("sample" %in% names(sample_table)) || !("R1" %in% names(sample_table))) {
    stop("\u274C samples.csv must contain at least 'sample' and 'R1' columns.")
  }
  
  has_R2 <- "R2" %in% names(sample_table)
  
  # Prepare output directory: results/merged/
  output_dir <- file.path(base_dir, "results", "merged")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Merge reads with PEAR
  commands <- character()
  for (i in seq_len(nrow(sample_table))) {
    sample <- sample_table$sample[i]
    r1 <- sample_table$R1[i]
    
    if (has_R2 && !is.na(sample_table$R2[i]) && sample_table$R2[i] != "") {
      r2 <- sample_table$R2[i]
      output_prefix <- file.path(output_dir, paste0(sample, "_ANC"))
      cmd <- sprintf("pear -f %s -r %s -o %s",
                     shQuote(r1), shQuote(r2), shQuote(output_prefix))
      message("\u25B6 Merging reads for sample: ", sample)
      status <- system(cmd)
      if (status != 0L) {
        stop("PEAR failed for sample ", sample,
             " with exit code ", status, ".")
      }
      commands <- c(commands, cmd)
    } else {
      message("\u2139 Skipping sample ", sample, " \u2014 only R1 found (single-end mode).")
    }
  }
  
  invisible(commands)
}
