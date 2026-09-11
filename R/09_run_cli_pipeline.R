#' Map Single-End or Overlapping Paired-End Reads
#'
#' Run FastQC, merge overlapping pairs with PEAR when R2 is supplied, and map
#' each sample with minimap2/samtools to a sorted, indexed BAM. R1-only samples
#' map directly. Extract barcodes from the returned BAMs with [barbac_xtr()].
#'
#' @param sample_table A data.frame with `sample`, `R1`, and optional `R2`, or
#'   the path to a CSV (or directory containing `samples.csv`). Missing, NA or
#'   empty R2 entries select single-end mode for that row. Sample names must
#'   be unique, start with a letter or digit, and contain only letters, digits,
#'   underscores, dots or hyphens. Read paths are relative to the working directory.
#' @param reference Path to the mapping-reference FASTA.
#' @param output_dir New or empty output directory. Existing results are never
#'   discovered as inputs or silently overwritten. Default: `"results"`.
#' @param verbose Print progress messages. Default: TRUE.
#' @param log_file New log-file path; defaults to `output_dir/pipeline.log`.
#' @param create_output_dir Create the output directory if needed. Default: TRUE.
#'
#' @details
#' PEAR is required only when at least one sample has R2. Paired mode maps the
#' assembled reads and excludes unmerged pairs. It is intended for overlapping
#' reads, not general paired-end mapping. A table may mix paired and R1-only rows.
#'
#' Mapping uses minimap2's short-read preset (`-x sr`) with secondary alignments
#' disabled. Secondary and supplementary alignments are removed from the BAM;
#' mapping statistics therefore count primary reads or merged molecules.
#' FastQC input filenames must produce unique report names within a run.
#'
#' Required commands stop the pipeline on failure, with their output recorded
#' in the log. MultiQC runs when available; its failure raises a warning and is
#' reported as `multiqc_status = "failed"`. This wrapper does not extract barcodes,
#' deduplicate UMIs, or apply study-specific filtering.
#'
#' @return An invisible list containing `commands`, `output_dir`, `fastqc_dir`,
#'   `merged_dir`, `bam_dir`, `stats`, `summary_file`, `log_file`, `multiqc_status`,
#'   and `samples`. The `samples` table links original sample labels to mode,
#'   mapping input and indexed BAM. `bam_files` is a vector named by sample.
#'   Paired BAM names retain `<sample>_ANC.assembled_sorted.bam`; R1-only BAMs
#'   use `<sample>_sorted.bam`. The `stats$sample` column retains these basenames
#'   without `_sorted.bam` for compatibility with earlier paired runs.
#' @md
#' @export
#' @examples
#' \dontrun{
#' samples <- data.frame(sample = "sample1", R1 = "sample1_R1.fastq.gz")
#' pipeline <- run_cli_pipeline(samples, "cassette.fasta", "results")
#' pipeline$bam_files[["sample1"]]
#' pipeline$stats
#' }
run_cli_pipeline <- function(sample_table, reference, output_dir = "results",
                             verbose = TRUE, log_file = NULL,
                             create_output_dir = TRUE) {
  scalar_path <- function(x, label) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
      stop(label, " must be a nonempty path.", call. = FALSE)
    path.expand(x)
  }
  if (is.character(sample_table)) {
    sample_path <- scalar_path(sample_table, "sample_table")
    if (dir.exists(sample_path)) sample_path <- file.path(sample_path, "samples.csv")
    if (!file.exists(sample_path)) stop("Sample table not found: ", sample_path, call. = FALSE)
    sample_table <- utils::read.csv(sample_path, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(sample_table) || !all(c("sample", "R1") %in% names(sample_table)) ||
      nrow(sample_table) == 0L || anyDuplicated(names(sample_table)))
    stop("sample_table must have at least one row and unique 'sample' and 'R1' columns.", call. = FALSE)
  labels <- as.character(sample_table$sample)
  if (anyNA(labels) || any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", labels)) ||
      anyDuplicated(tolower(labels)))
    stop("Sample names must be unique (ignoring case), start with a letter or digit, and contain only letters, digits, '.', '_' or '-'.", call. = FALSE)
  r1 <- as.character(sample_table$R1)
  r2 <- if ("R2" %in% names(sample_table)) as.character(sample_table$R2) else rep(NA_character_, length(r1))
  paired <- !is.na(r2) & nzchar(trimws(r2))
  read_paths <- c(r1, r2[paired])
  if (anyNA(read_paths) || any(!nzchar(read_paths)) ||
      any(!file.exists(read_paths) | dir.exists(read_paths)))
    stop("Every supplied R1/R2 must point to an existing FASTQ file.", call. = FALSE)
  r1 <- normalizePath(r1, mustWork = TRUE)
  r2[paired] <- normalizePath(r2[paired], mustWork = TRUE)
  read_paths <- unique(c(r1, r2[paired]))
  qc_names <- tolower(sub("[.](fastq|fq)([.]gz)?$", "", basename(read_paths), ignore.case = TRUE))
  if (anyDuplicated(qc_names))
    stop("FASTQ filenames must produce unique FastQC report names; rename files with repeated basenames.", call. = FALSE)
  reference <- scalar_path(reference, "reference")
  if (!file.exists(reference) || dir.exists(reference)) stop("Reference file not found: ", reference, call. = FALSE)
  reference <- normalizePath(reference, mustWork = TRUE)
  bases <- ifelse(paired, paste0(labels, "_ANC.assembled"), labels)
  if (anyDuplicated(tolower(bases))) stop("Sample names produce colliding BAM filenames.", call. = FALSE)
  output_dir <- scalar_path(output_dir, "output_dir")
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    stop("output_dir must be new or empty; choose a new directory to preserve existing results.", call. = FALSE)
  if (!dir.exists(output_dir) && !isTRUE(create_output_dir))
    stop("Output directory does not exist; set create_output_dir = TRUE.", call. = FALSE)
  required <- c("fastqc", if (any(paired)) "pear", "minimap2", "samtools")
  bins <- .barbac_require_tools(required)
  multiqc_bin <- .barbac_tool("multiqc")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  if (is.null(log_file)) log_file <- file.path(output_dir, "pipeline.log")
  log_file <- scalar_path(log_file, "log_file")
  if (file.exists(log_file)) stop("Log file already exists: ", log_file, call. = FALSE)
  log_conn <- file(log_file, open = "wt")
  on.exit(close(log_conn), add = TRUE)
  log_msg <- function(msg) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = log_conn)
    flush(log_conn)
    if (verbose) message(msg)
  }
  commands <- character()
  run <- function(tool, args, stdout = NULL, required = TRUE) {
    bin <- if (tool == "multiqc") multiqc_bin else bins[[tool]]
    cmd <- paste(c(shQuote(bin), shQuote(as.character(args))), collapse = " ")
    if (!is.null(stdout)) cmd <- paste(cmd, ">", shQuote(stdout))
    commands <<- c(commands, cmd)
    log_msg(cmd)
    diagnostic <- tempfile("barbac-command-")
    on.exit(unlink(diagnostic), add = TRUE)
    status <- system2(bin, shQuote(as.character(args)),
      stdout = if (is.null(stdout)) diagnostic else stdout, stderr = diagnostic)
    if (file.exists(diagnostic)) {
      cat(readLines(diagnostic, warn = FALSE), sep = "\n", file = log_conn)
      flush(log_conn)
    }
    if (status != 0L) {
      msg <- paste(tool, "failed with exit code", status, "- see", log_file)
      log_msg(msg)
      if (required) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
    }
    status
  }
  log_msg(paste("BARBAC PIPELINE:", length(labels), "samples;", sum(paired), "paired;", sum(!paired), "R1-only"))
  fastqc_dir <- file.path(output_dir, "fastQC")
  merged_dir <- file.path(output_dir, "merged")
  bam_dir <- file.path(merged_dir, "bam")
  for (path in c(fastqc_dir, merged_dir, bam_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  for (fq in read_paths) run("fastqc", c(fq, "-o", fastqc_dir))
  mapping_files <- r1
  for (i in which(paired)) {
    prefix <- file.path(merged_dir, paste0(labels[i], "_ANC"))
    run("pear", c("-f", r1[i], "-r", r2[i], "-o", prefix))
    mapping_files[i] <- paste0(prefix, ".assembled.fastq")
    if (!file.exists(mapping_files[i]) || file.info(mapping_files[i])$size == 0)
      stop("PEAR produced no assembled reads for sample ", labels[i], "; overlapping pairs are required.", call. = FALSE)
  }
  bam_files <- stats::setNames(file.path(bam_dir, paste0(bases, "_sorted.bam")), labels)
  stats <- vector("list", length(labels))
  for (i in seq_along(labels)) {
    sam <- file.path(bam_dir, paste0(bases[i], ".sam"))
    primary <- file.path(bam_dir, paste0(bases[i], ".bam"))
    run("minimap2", c("-a", "-x", "sr", "--secondary=no", reference, mapping_files[i]), stdout = sam)
    run("samtools", c("view", "-b", "-F", "2304", "-o", primary, sam))
    run("samtools", c("sort", "-o", bam_files[i], primary))
    run("samtools", c("index", bam_files[i]))
    if (!file.exists(bam_files[i]) || !file.exists(paste0(bam_files[i], ".bai")))
      stop("Missing indexed BAM for sample ", labels[i], call. = FALSE)
    count_file <- tempfile("barbac-count-")
    on.exit(unlink(count_file), add = TRUE)
    read_count <- function(args) {
      run("samtools", c("view", "-c", args, bam_files[i]), stdout = count_file)
      value <- suppressWarnings(as.numeric(readLines(count_file, warn = FALSE)))
      if (length(value) != 1L || !is.finite(value) || value < 0 || value != floor(value))
        stop("Invalid samtools count for sample ", labels[i], call. = FALSE)
      value
    }
    stats[[i]] <- tibble::tibble(sample = bases[i], mapped = read_count(c("-F", "2308")),
                                unmapped = read_count(c("-f", "4", "-F", "2304")))
    unlink(c(sam, primary, count_file))
  }
  stats_df <- dplyr::bind_rows(stats)
  summary_file <- file.path(output_dir, "bam_summary.csv")
  readr::write_csv(stats_df, summary_file)
  multiqc_status <- "unavailable"
  if (nzchar(multiqc_bin)) {
    mqc_dir <- file.path(output_dir, "multiqc")
    dir.create(mqc_dir, showWarnings = FALSE)
    status <- run("multiqc", c(fastqc_dir, "-o", mqc_dir), required = FALSE)
    multiqc_status <- if (status == 0L) "completed" else "failed"
  }
  samples <- tibble::tibble(sample = labels, mode = ifelse(paired, "paired", "single"),
    mapping_input = mapping_files, bam_file = unname(bam_files))
  readr::write_csv(samples, file.path(output_dir, "sample_outputs.csv"))
  log_msg(paste("Mapping complete:", sum(stats_df$mapped), "mapped;", sum(stats_df$unmapped), "unmapped."))
  invisible(list(commands = commands, output_dir = output_dir, fastqc_dir = fastqc_dir,
    merged_dir = merged_dir, bam_dir = bam_dir, stats = stats_df, summary_file = summary_file,
    log_file = normalizePath(log_file), samples = samples, bam_files = bam_files,
    multiqc_status = multiqc_status))
}
