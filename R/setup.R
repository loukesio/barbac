#' Configure Conda Environment for barbac
#'
#' This function creates and configures a conda environment for the barbac pipeline,
#' installing necessary bioinformatics command-line tools like fastqc, multiqc, PEAR, minimap2, and samtools.
#'
#' @param env_name Name of the conda environment to create. Default is "barbac_env".
#' @param channels List of conda channels to use (default: c("conda-forge", "bioconda"))
#' @param confirm If FALSE, it will skip the interactive menu and install directly. Default is TRUE.
#'
#' @return Prints progress messages and returns the environment name if successful.
#' @export
configure_environment <- function(env_name = "barbac_env",
                                  channels = c("conda-forge", "bioconda"),
                                  confirm = TRUE) {
  
  # Check for conda
  if (is.null(reticulate::conda_version())) {
    stop("Conda is not installed on your system. Please install Conda first.")
  }
  
  # Check existing envs
  envs <- reticulate::conda_list()
  if (env_name %in% envs$name) {
    message(sprintf("The conda environment '%s' already exists. Skipping creation.", env_name))
    return(invisible(env_name))
  }
  
  if (interactive() && confirm) {
    choice <- utils::menu(c("Yes", "No"),
                          title = sprintf("The conda environment '%s' does not exist. Create it?", env_name))
    if (choice != 1) stop("Environment creation cancelled by user.")
  }
  
  # Required CLI tools
  tools <- c("fastqc", "multiqc", "pear", "minimap2", "samtools")
  
  message("Creating conda environment and installing tools: ", paste(tools, collapse = ", "))
  
  reticulate::conda_create(
    envname = env_name,
    packages = tools,
    channel = channels
  )
  
  message("✅ Environment '", env_name, "' created successfully.")
  message("🎉 You can now use barbac's functions like run_fastqc(), run_pear_merge(), etc.")
  
  invisible(env_name)
}
