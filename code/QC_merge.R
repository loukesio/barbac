#' Run qc_merge
#'
#' This function runs fastqc and bbmerge on the given data directory based on the sample table.
#' If the tools are not available, it runs the commands using a Docker container.
#'
#' @param sample_table Path to the CSV file containing sample information
#' @param data_dir Directory containing the data files
#' @param config_file Optional configuration file in YAML format
#' @export
qc_merge <- function(sample_table, data_dir, config_file = NULL) {
  # Load necessary libraries
  library(sys)
  library(readr)
  library(yaml)
  
  # Read the sample table
  samples <- read_csv(sample_table)
  
  # Check if necessary columns are present
  required_columns <- c("sample_names", "R1_files", "R2_files")
  if (!all(required_columns %in% colnames(samples))) {
    stop("Sample table must contain the columns: sample_names, R1_files, and R2_files.")
  }
  
  # Define output directories
  fastqc_output_dir <- file.path(data_dir, "fastqc_output")
  merged_dir <- file.path(data_dir, "merged")
  
  # Create output directories if they don't exist
  if (!dir.exists(fastqc_output_dir)) {
    dir.create(fastqc_output_dir, recursive = TRUE)
  }
  if (!dir.exists(merged_dir)) {
    dir.create(merged_dir, recursive = TRUE)
  }
  
  # Read configuration file if provided
  config <- list()
  if (!is.null(config_file)) {
    config <- yaml::read_yaml(config_file)
  }
  
  # Function to check if a command exists
  command_exists <- function(command) {
    sys::exec_internal("which", command)$status == 0
  }
  
  # Function to run Docker container
  run_docker <- function() {
    docker_command <- sprintf(
      "docker run -v %s:/data your_dockerhub_username/barbac_pipeline:v1.2 Rscript -e 'barbac::qc_merge(\"/data/%s\", \"/data\", \"%s\")'",
      normalizePath(data_dir), basename(sample_table), ifelse(is.null(config_file), "", basename(config_file))
    )
    message("Running Docker command: ", docker_command)
    system(docker_command)
  }
  
  # Check if fastqc and bbmerge are installed
  if (!command_exists("fastqc") || !command_exists("bbmerge.sh")) {
    message("fastqc or bbmerge.sh not found, running in Docker container...")
    run_docker()
    return()
  }
  
  # Run fastqc on each R1 and R2 file
  fastqc_threads <- ifelse(!is.null(config$fastqc$threads), config$fastqc$threads, 1)
  
  for (i in 1:nrow(samples)) {
    r1_file <- file.path(data_dir, samples$R1_files[i])
    r2_file <- file.path(data_dir, samples$R2_files[i])
    sys::exec_wait("fastqc", c(r1_file, "-o", fastqc_output_dir, "--threads", fastqc_threads))
    sys::exec_wait("fastqc", c(r2_file, "-o", fastqc_output_dir, "--threads", fastqc_threads))
  }
  
  # Run bbmerge.sh on each pair of fastq files
  for (i in 1:nrow(samples)) {
    sample <- samples$sample_names[i]
    r1 <- file.path(data_dir, samples$R1_files[i])
    r2 <- file.path(data_dir, samples$R2_files[i])
    merged_file <- file.path(merged_dir, paste0(sample, "_merged.fastq"))
    unmerged1_file <- file.path(merged_dir, paste0(sample, "_unmerged1.fastq"))
    unmerged2_file <- file.path(merged_dir, paste0(sample, "_unmerged2.fastq"))
    
    bbmerge_params <- list(
      in1 = r1,
      in2 = r2,
      out = merged_file,
      outu1 = unmerged1_file,
      outu2 = unmerged2_file
    )
    
    sys::exec_wait("bbmerge.sh", unlist(bbmerge_params))
  }
  
  # Handle config file if provided
  if (!is.null(config_file)) {
    file.copy(config_file, merged_dir)
  }
}
