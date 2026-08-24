#' Run Full Barbac CLI Pipeline
#'
#' This function runs the full CLI pipeline for barcode analysis: FastQC -> PEAR -> Minimap2 -> BAM Stats.
#'
#' @param sample_table Path to samples.csv or a data.frame/tibble containing `sample`, `R1`, and optionally `R2`.
#' @param reference Path to reference FASTA file.
#' @param output_dir Directory to write output files. Default: "results".
#' @param verbose Whether to print progress messages to console. Default: TRUE.
#' @param log_file File path for logging. If NULL, uses output_dir/pipeline.log. Default: NULL.
#' @param create_output_dir Whether to create output_dir if it doesn't exist. Default: TRUE.
#'
#' @return Invisible list containing:
#'   \itemize{
#'     \item commands: Character vector of system commands executed
#'     \item output_dir: Path to output directory
#'     \item stats: Data frame with BAM statistics
#'   }
#' @importFrom utils read.csv txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' \dontrun{
#' # Run with default output directory (./results)
#' run_cli_pipeline(
#'   sample_table = "samples.csv",
#'   reference = "barcodes.fasta"
#' )
#'
#' # Run with custom output directory
#' run_cli_pipeline(
#'   sample_table = "samples.csv",
#'   reference = "barcodes.fasta",
#'   output_dir = "/path/to/my_results"
#' )
#' }
run_cli_pipeline <- function(sample_table,
                             reference,
                             output_dir = "results",
                             verbose = TRUE,
                             log_file = NULL,
                             create_output_dir = TRUE) {
  
  # -----------------------------
  # Validate and create output directory
  # -----------------------------
  if (create_output_dir) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      if (verbose) message("\u2713 Created output directory: ", output_dir)
    }
  } else {
    if (!dir.exists(output_dir)) {
      stop("\u274C Output directory does not exist: ", output_dir, 
           "\n   Set create_output_dir = TRUE to create it automatically.")
    }
  }
  
  # Make output_dir absolute path for clarity
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  
  # Set default log file location
  if (is.null(log_file)) {
    log_file <- file.path(output_dir, "pipeline.log")
  }
  
  # -----------------------------
  # Handle sample_table
  # -----------------------------
  if (is.character(sample_table)) {
    if (dir.exists(sample_table)) {
      sample_path <- file.path(sample_table, "samples.csv")
    } else {
      sample_path <- sample_table
    }
    if (!file.exists(sample_path)) {
      stop("\u274C Could not find samples.csv at: ", sample_path)
    }
    sample_table <- read.csv(sample_path)
  }
  
  if (!("sample" %in% names(sample_table)) || !("R1" %in% names(sample_table))) {
    stop("\u274C samples.csv must have at least 'sample' and 'R1' columns.")
  }
  
  # -----------------------------
  # Validate reference file
  # -----------------------------
  if (!file.exists(reference)) {
    stop("\u274C Reference file not found: ", reference)
  }
  
  # -----------------------------
  # Prepare logging
  # -----------------------------
  log_conn <- NULL
  try({
    log_conn <- file(log_file, open = "wt")
  }, silent = TRUE)
  
  if (is.null(log_conn) || !isOpen(log_conn)) {
    stop("\u274C Failed to open log file: ", log_file)
  }
  
  on.exit({
    if (!is.null(log_conn) && isOpen(log_conn)) close(log_conn)
  }, add = TRUE)
  
  log_msg <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- sprintf("[%s] %s\n", timestamp, msg)
    try(cat(line, file = log_conn, append = TRUE), silent = TRUE)
    if (verbose) message(msg)
  }
  
  # Log initial setup
  log_msg("=" %R% 60)
  log_msg("BARBAC PIPELINE STARTED")
  log_msg("=" %R% 60)
  log_msg(paste("Output directory:", output_dir))
  log_msg(paste("Reference file:", reference))
  log_msg(paste("Number of samples:", nrow(sample_table)))
  log_msg("")
  
  # -----------------------------
  # Start progress
  # -----------------------------
  steps <- c("FastQC", "PEAR", "minimap2", "BAM stats")
  if (verbose) pb <- txtProgressBar(min = 0, max = length(steps), style = 3)
  all_cmds <- character()
  
  # -----------------------------
  # Step 1: FastQC
  # -----------------------------
  log_msg("\u25B6 STEP 1/4: Running FastQC")
  fastqc_dir <- file.path(output_dir, "fastQC")
  dir.create(fastqc_dir, recursive = TRUE, showWarnings = FALSE)
  
  fastq_files <- sample_table$R1
  if ("R2" %in% names(sample_table)) {
    fastq_files <- unique(c(fastq_files, sample_table$R2[!is.na(sample_table$R2)]))
  }
  
  fastqc_cmds <- vapply(fastq_files, function(fq) {
    if (!file.exists(fq)) {
      log_msg(paste("\u26A0 Warning: FASTQ file not found:", fq))
      return("")
    }
    cmd <- sprintf("fastqc %s -o %s", shQuote(fq), shQuote(fastqc_dir))
    log_msg(paste("  Running:", cmd))
    exit_code <- system(cmd)
    if (exit_code != 0) {
      log_msg(paste("  \u26A0 FastQC failed with exit code:", exit_code))
    }
    cmd
  }, character(1))
  
  all_cmds <- c(all_cmds, fastqc_cmds[fastqc_cmds != ""])
  if (verbose) setTxtProgressBar(pb, 1)
  log_msg(paste("\u2713 FastQC completed. Results in:", fastqc_dir))
  log_msg("")
  
  # -----------------------------
  # Step 2: PEAR
  # -----------------------------
  log_msg("\u25B6 STEP 2/4: Merging reads with PEAR")
  merged_dir <- file.path(output_dir, "merged")
  dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)
  
  pear_cmds <- character()
  for (i in seq_len(nrow(sample_table))) {
    sample <- sample_table$sample[i]
    r1 <- sample_table$R1[i]
    
    if (!file.exists(r1)) {
      log_msg(paste("  \u26A0 Skipping sample", sample, "- R1 file not found:", r1))
      next
    }
    
    if ("R2" %in% names(sample_table) && 
        !is.na(sample_table$R2[i]) && 
        sample_table$R2[i] != "") {
      r2 <- sample_table$R2[i]
      
      if (!file.exists(r2)) {
        log_msg(paste("  \u26A0 Skipping sample", sample, "- R2 file not found:", r2))
        next
      }
      
      out_pref <- file.path(merged_dir, paste0(sample, "_ANC"))
      cmd <- sprintf("pear -f %s -r %s -o %s", 
                     shQuote(r1), shQuote(r2), shQuote(out_pref))
      log_msg(paste("  Merging sample:", sample))
      exit_code <- system(cmd)
      if (exit_code != 0) {
        log_msg(paste("  \u26A0 PEAR failed with exit code:", exit_code))
      }
      pear_cmds <- c(pear_cmds, cmd)
    } else {
      log_msg(paste("  \u26A0 Skipping PEAR merge for sample", sample, "(single-end or R2 missing)"))
    }
  }
  
  all_cmds <- c(all_cmds, pear_cmds)
  if (verbose) setTxtProgressBar(pb, 2)
  log_msg(paste("\u2713 PEAR completed. Results in:", merged_dir))
  log_msg("")
  
  # -----------------------------
  # Step 3: Minimap2
  # -----------------------------
  log_msg("\u25B6 STEP 3/4: Mapping with minimap2 + samtools")
  bam_dir <- file.path(merged_dir, "bam")
  dir.create(bam_dir, recursive = TRUE, showWarnings = FALSE)
  
  mapping_inputs <- .pipeline_mapping_inputs(sample_table, merged_dir)
  merged_files <- mapping_inputs$merged_files
  single_files <- mapping_inputs$single_files
  mapping_files <- mapping_inputs$files
  mapping_names <- mapping_inputs$names

  if (length(mapping_files) == 0) {
    stop("\u274C No merged paired-end or single-end FASTQ files available for mapping.")
  }
  
  log_msg(paste("  Found", length(mapping_files), "FASTQ files to process (",
                length(merged_files), "merged,", length(single_files), "single-end)"))
  
  map_cmds <- vapply(seq_along(mapping_files), function(i) {
    fq <- mapping_files[i]
    base <- mapping_names[i]
    sam <- file.path(bam_dir, paste0(base, ".sam"))
    bam <- file.path(bam_dir, paste0(base, ".bam"))
    sorted_bam <- file.path(bam_dir, paste0(base, "_sorted.bam"))
    
    cmds <- c(
      sprintf("minimap2 -a %s %s > %s", 
              shQuote(reference), shQuote(fq), shQuote(sam)),
      sprintf("samtools view -S -b %s > %s", 
              shQuote(sam), shQuote(bam)),
      sprintf("samtools sort %s -o %s", 
              shQuote(bam), shQuote(sorted_bam)),
      sprintf("samtools index %s", 
              shQuote(sorted_bam))
    )
    
    log_msg(paste("  Processing:", basename(fq)))
    for (cmd in cmds) {
      exit_code <- system(cmd)
      if (exit_code != 0) {
        log_msg(paste("  \u26A0 Command failed with exit code:", exit_code))
      }
    }
    
    # Clean up intermediate files
    if (file.exists(sam)) file.remove(sam)
    if (file.exists(bam)) file.remove(bam)
    
    paste(cmds, collapse = " && ")
  }, character(1))
  
  all_cmds <- c(all_cmds, map_cmds)
  if (verbose) setTxtProgressBar(pb, 3)
  log_msg(paste("\u2713 Minimap2 completed. Results in:", bam_dir))
  log_msg("")
  
  # -----------------------------
  # Step 4: BAM Stats
  # -----------------------------
  log_msg("\u25B6 STEP 4/4: Summarising BAM stats")
  bam_files <- list.files(bam_dir, pattern = "_sorted\\.bam$", full.names = TRUE)
  
  if (length(bam_files) == 0) {
    warning("No sorted BAM files found for statistics")
    stats_df <- data.frame(sample = character(), mapped = integer(), unmapped = integer())
  } else {
    stats <- lapply(bam_files, function(bam) {
      sample <- gsub("_sorted\\.bam$", "", basename(bam))
      mapped <- as.integer(system2("samtools", c("view", "-c", "-F", "4", bam), stdout = TRUE))
      unmapped <- as.integer(system2("samtools", c("view", "-c", "-f", "4", bam), stdout = TRUE))
      tibble::tibble(sample, mapped, unmapped)
    })
    stats_df <- dplyr::bind_rows(stats)
  }
  
  summary_file <- file.path(output_dir, "bam_summary.csv")
  readr::write_csv(stats_df, summary_file)
  log_msg(paste("\u2713 BAM stats saved to:", summary_file))
  
  if (verbose) setTxtProgressBar(pb, 4)
  log_msg("")
  
  # -----------------------------
  # Step 5: MultiQC (optional)
  # -----------------------------
  if (Sys.which("multiqc") != "") {
    log_msg("\u25B6 Running MultiQC on FastQC results...")
    mqc_dir <- file.path(output_dir, "multiqc")
    dir.create(mqc_dir, recursive = TRUE, showWarnings = FALSE)
    mqc_cmd <- sprintf("multiqc %s -o %s", shQuote(fastqc_dir), shQuote(mqc_dir))
    system(mqc_cmd)
    log_msg(paste("\u2713 MultiQC report generated in:", mqc_dir))
    all_cmds <- c(all_cmds, mqc_cmd)
  } else {
    log_msg("\u26A0 MultiQC not found in PATH \u2014 skipping.")
  }
  
  # -----------------------------
  # Summary
  # -----------------------------
  if (verbose) close(pb)
  log_msg("")
  log_msg("=" %R% 60)
  log_msg("PIPELINE SUMMARY")
  log_msg("=" %R% 60)
  log_msg(paste("Total samples processed:", nrow(sample_table)))
  log_msg(paste("Merged files created:", length(merged_files)))
  log_msg(paste("BAM files created:", length(bam_files)))
  if (nrow(stats_df) > 0) {
    log_msg(paste("Total mapped reads:", sum(stats_df$mapped)))
    log_msg(paste("Total unmapped reads:", sum(stats_df$unmapped)))
  }
  log_msg("")
  log_msg("OUTPUT STRUCTURE:")
  log_msg(paste("\u251C\u2500\u2500 FastQC results:", fastqc_dir))
  log_msg(paste("\u251C\u2500\u2500 Merged reads:", merged_dir))
  log_msg(paste("\u251C\u2500\u2500 BAM files:", bam_dir))
  log_msg(paste("\u2514\u2500\u2500 Summary:", summary_file))
  log_msg("")
  log_msg("\u2705 Barbac pipeline completed successfully!")
  log_msg("=" %R% 60)
  
  # Return results invisibly
  invisible(list(
    commands = all_cmds,
    output_dir = output_dir,
    fastqc_dir = fastqc_dir,
    merged_dir = merged_dir,
    bam_dir = bam_dir,
    stats = stats_df,
    summary_file = summary_file,
    log_file = log_file
  ))
}

# Helper operator for string repetition (for log formatting)
`%R%` <- function(x, n) paste(rep(x, n), collapse = "")

# Internal: resolve FASTQ inputs for the mapping stage. Paired-end samples are
# represented by PEAR's assembled output; samples without R2 go straight from
# R1 to minimap2.
.pipeline_mapping_inputs <- function(sample_table, merged_dir) {
  merged_files <- list.files(merged_dir,
                             pattern = "([_.])assembled\\.fastq$",
                             full.names = TRUE)
  has_r2 <- if ("R2" %in% names(sample_table)) {
    !is.na(sample_table$R2) & sample_table$R2 != ""
  } else {
    rep(FALSE, nrow(sample_table))
  }
  single_rows <- which(!has_r2 & file.exists(sample_table$R1))
  single_files <- sample_table$R1[single_rows]

  list(
    files = c(merged_files, single_files),
    names = c(tools::file_path_sans_ext(basename(merged_files)),
              as.character(sample_table$sample[single_rows])),
    merged_files = merged_files,
    single_files = single_files
  )
}
