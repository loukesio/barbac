#' Check availability of tools in barbac conda environment
#'
#' @param env_name Name of the conda environment (default: "barbac_env")
#' @param tools Character vector of tool names to check
#'
#' @return A data.frame with tool names and availability
#' @export
check_barbac_tools <- function(env_name = "barbac_env",
                               tools = c("fastqc", "multiqc", "pear", "minimap2", "samtools")) {
  
  results <- lapply(tools, function(tool) {
    cmd <- sprintf("bash -c 'source activate %s && which %s 2>/dev/null'", env_name, tool)
    path <- suppressWarnings(system(cmd, intern = TRUE, ignore.stderr = TRUE))
    
    available <- length(path) > 0 && nzchar(path[1])
    
    version <- NA_character_
    if (available) {
      version_cmd <- sprintf("bash -c 'source activate %s && %s --version 2>&1 | head -n1'", 
                             env_name, tool)
      version <- tryCatch(
        system(version_cmd, intern = TRUE)[1],
        error = function(e) NA_character_
      )
    }
    
    data.frame(
      tool = tool,
      available = available,
      version = ifelse(is.na(version), NA_character_, version),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, results)
}