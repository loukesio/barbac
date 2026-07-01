#' Configure Conda Environment for barbac
#'
#' @param env_name Name of the conda environment to create. Default is "barbac_env".
#' @param channels List of conda channels to use (default: c("conda-forge", "bioconda"))
#' @param confirm If FALSE, it will skip the interactive menu. Default is TRUE.
#'
#' @return Prints progress messages and returns the environment path if successful.
#' @export
configure_environment <- function(env_name = "barbac_env",
                                  channels = c("conda-forge", "bioconda"),
                                  confirm = TRUE) {
  
  conda_exe <- reticulate::conda_binary()
  if (is.null(conda_exe)) {
    stop("Conda is not installed. Please install Conda first.")
  }
  
  message("Using conda: ", conda_exe)
  
  # Check existing environments
  env_list <- system2(conda_exe, c("env", "list"), stdout = TRUE, stderr = TRUE)
  env_exists <- any(grepl(paste0("\\s", env_name, "\\s+"), env_list) | 
                      grepl(paste0("^", env_name, "\\s+"), env_list))
  
  if (env_exists) {
    env_line <- grep(env_name, env_list, value = TRUE)[1]
    env_path <- sub(".*\\s(/.*)", "\\1", env_line)
    
    message("\u2705 Environment '", env_name, "' already exists at:")
    message("   ", env_path)
    message("\U0001F389 Ready to use barbac functions.")
    return(invisible(env_path))
  }
  
  if (interactive() && confirm) {
    choice <- utils::menu(c("Yes", "No"),
                          title = sprintf("Create conda environment '%s'?", env_name))
    if (choice != 1) {
      message("Cancelled.")
      return(invisible(NULL))
    }
  }
  
  tools <- c("fastqc", "multiqc", "pear", "minimap2", "samtools")
  
  message("Creating environment with tools: ", paste(tools, collapse = ", "))
  message("This may take several minutes...")
  
  channel_args <- unlist(lapply(channels, function(ch) c("-c", ch)))
  
  result <- system2(
    conda_exe,
    c("create", "--yes", "--name", env_name, tools, channel_args),
    stdout = TRUE,
    stderr = TRUE
  )
  
  Sys.sleep(2)
  env_list_after <- system2(conda_exe, c("env", "list"), stdout = TRUE, stderr = TRUE)
  env_created <- any(grepl(paste0("\\s", env_name, "\\s+"), env_list_after))
  
  if (env_created) {
    env_line <- grep(env_name, env_list_after, value = TRUE)[1]
    env_path <- sub(".*\\s(/.*)", "\\1", env_line)
    message("\u2705 Environment created at: ", env_path)
    return(invisible(env_path))
  } else {
    warning("Creation may have failed:\n", paste(result, collapse = "\n"))
    return(invisible(NULL))
  }
}