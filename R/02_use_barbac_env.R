#' Use barbac conda environment in the current R session
#'
#' Sets PATH to include the barbac conda environment so CLI tools work in R.
#'
#' @param env_name Name of the conda environment. Default is "barbac_env"
#'
#' @return The full path to the environment's bin directory, invisibly
#' @export
use_barbac_env <- function(env_name = "barbac_env") {
  
  conda_exe <- reticulate::conda_binary()
  if (is.null(conda_exe)) {
    stop("Conda not found. Please install Conda first.")
  }
  
  # Get list of environments
  env_list <- system2(conda_exe, c("env", "list"), stdout = TRUE, stderr = TRUE)
  
  # Find lines containing the environment name
  env_line <- grep(paste0("\\b", env_name, "\\b"), env_list, value = TRUE)
  
  if (length(env_line) == 0) {
    stop(sprintf("Environment '%s' not found. Run configure_environment() first.", env_name))
  }
  
  # Parse the path: split by whitespace and find paths starting with /
  parts <- strsplit(env_line[1], "\\s+")[[1]]
  parts <- parts[nchar(parts) > 0]  # Remove empty strings
  
  # Find the path (starts with /)
  path_candidates <- parts[grepl("^/", parts)]
  
  if (length(path_candidates) == 0) {
    stop(sprintf("Could not parse environment path from: %s", env_line[1]))
  }
  
  env_path <- path_candidates[length(path_candidates)]
  
  # Remove any trailing markers like * or +
  env_path <- gsub("[*+]\\s*$", "", env_path)
  env_path <- trimws(env_path)
  
  # Verify the environment path exists
  if (!dir.exists(env_path)) {
    stop(sprintf("Environment directory does not exist: '%s'", env_path))
  }
  
  # Construct bin path
  bin_path <- file.path(env_path, "bin")
  
  if (!dir.exists(bin_path)) {
    stop(sprintf("Bin directory not found at: '%s'\nThe environment may not be properly configured.", bin_path))
  }
  
  # Update PATH environment variable
  current_path <- Sys.getenv("PATH")
  
  if (!grepl(bin_path, current_path, fixed = TRUE)) {
    new_path <- paste(bin_path, current_path, sep = .Platform$path.sep)
    Sys.setenv(PATH = new_path)
    message("\u2705 PATH updated to include: ", bin_path)
    message("   Bioinformatics tools are now available in this R session")
    message("   Available tools: fastqc, multiqc, pear, minimap2, samtools")
  } else {
    message("\U0001F501 PATH already includes: ", bin_path)
    message("   Bioinformatics tools are already available")
  }
  
  invisible(bin_path)
}
