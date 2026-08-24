#' Check availability of tools in barbac conda environment
#'
#' @param env_name Name of the conda environment (default: "barbac_env")
#' @param tools Character vector of tool names to check
#'
#' @return A data.frame with tool names and availability
#' @export
check_barbac_tools <- function(env_name = "barbac_env",
                               tools = c("fastqc", "multiqc", "pear", "minimap2", "samtools")) {
  conda_exe <- reticulate::conda_binary()
  if (is.null(conda_exe) || !nzchar(conda_exe)) {
    stop("Conda is not installed or could not be located.")
  }

  results <- lapply(tools, function(tool) {
    output <- tryCatch(
      suppressWarnings(system2(
        conda_exe,
        c("run", "--name", shQuote(env_name), shQuote(tool), "--version"),
        stdout = TRUE, stderr = TRUE
      )),
      error = function(e) structure(character(), status = 1L)
    )
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    available <- identical(as.integer(status), 0L)
    version <- if (available && length(output)) output[[1L]] else NA_character_
    
    data.frame(
      tool = tool,
      available = available,
      version = ifelse(is.na(version), NA_character_, version),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, results)
}
