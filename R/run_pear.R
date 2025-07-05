#' Merge paired-end reads using PEAR
#'
#' @param sample_table A data.frame or tibble with `sample`, `R1`, and `R2` columns
#' @param output_dir Where merged reads will be stored. Default is "merged"
#'
#' @return Character vector of PEAR commands run
#' @export
run_pear_merge <- function(sample_table, output_dir = "merged") {
  if (!all(c("sample", "R1", "R2") %in% colnames(sample_table))) {
    stop("sample_table must contain 'sample', 'R1', and 'R2' columns.")
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  commands <- vapply(seq_len(nrow(sample_table)), function(i) {
    sample <- sample_table$sample[i]
    r1 <- sample_table$R1[i]
    r2 <- sample_table$R2[i]
    output_prefix <- file.path(output_dir, paste0(sample, "_ANC"))
    cmd <- sprintf("pear -f %s -r %s -o %s", r1, r2, output_prefix)
    message("Running: ", cmd)
    system(cmd)
    return(cmd)
  }, character(1))
  
  invisible(commands)
}
