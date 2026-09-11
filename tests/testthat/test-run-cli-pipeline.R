cli_read_fixture <- function() {
  folder <- tempfile("barbac cli ' fixture ")
  dir.create(folder)
  set.seed(20260911)
  dna <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  left <- dna(170); right <- dna(170)
  parents <- replicate(3, dna(26))
  variants <- c(paste0(substr(parents[1], 1, 12), "A", substr(parents[1], 13, 26)),
                paste0(substr(parents[1], 1, 12), substr(parents[1], 14, 26)))
  truth <- data.frame(barcode = c(parents, variants), counts = c(80, 40, 20, 2, 2))
  seqs <- paste0(left, rep(truth$barcode, truth$counts), right)
  write_reads <- function(path, sequences, mate) {
    lines <- unlist(lapply(seq_along(sequences), function(i)
      c(paste0("@read", i, "/", mate), sequences[i], "+", strrep("I", nchar(sequences[i])))))
    con <- gzfile(path, "wt"); on.exit(close(con)); writeLines(lines, con)
  }
  r1 <- file.path(folder, "reads_R1.fastq.gz"); r2 <- file.path(folder, "reads_R2.fastq.gz")
  write_reads(r1, substr(seqs, 1, 250), 1)
  write_reads(r2, substr(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs))), 1, 250), 2)
  reference <- file.path(folder, "reference.fasta")
  writeLines(c(">cassette", paste0(left, strrep("N", 26), right)), reference)
  list(folder = folder, r1 = r1, r2 = r2, reference = reference, truth = truth,
       pattern = paste0(substr(left, 159, 170), "([ACGT]{24,28})", substr(right, 1, 12)))
}

cli_tools_available <- function(paired = FALSE) {
  names <- c("fastqc", "minimap2", "samtools", if (paired) "pear")
  skip_if(any(!nzchar(vapply(names, .barbac_tool, character(1)))), "CLI executables unavailable")
}

cli_expect_extracted <- function(pipe, fixture) {
  for (bam in pipe$bam_files) {
    output <- barbac_xtr(bam, "cassette", 171, 196, flank_pattern = fixture$pattern, verbose = FALSE)
    counts <- read.csv(output)
    expect_setequal(counts$barcode, fixture$truth$barcode)
    expect_equal(counts$counts[match(fixture$truth$barcode, counts$barcode)], fixture$truth$counts)
    clustered <- super_cluster2(counts, verbose = FALSE)
    expect_equal(nrow(clustered), 3)
    expect_equal(sum(clustered$sum_counts), 144)
  }
}

test_that("R1-only reads map and extract exactly without requiring PEAR", {
  cli_tools_available()
  f <- cli_read_fixture(); on.exit(unlink(f$folder, recursive = TRUE))
  resolver <- .barbac_require_tools; lookup <- .barbac_tool
  local_mocked_bindings(
    .barbac_require_tools = function(tools, ...) { expect_false("pear" %in% tools); resolver(tools, ...) },
    .barbac_tool = function(tool, ...) if (tool == "multiqc") "" else lookup(tool, ...))
  pipe <- run_cli_pipeline(data.frame(sample = "single", R1 = f$r1), f$reference,
                           file.path(f$folder, "single results"), verbose = FALSE)
  expect_identical(pipe$samples$mode, "single")
  expect_equal(pipe$stats$mapped, 144)
  expect_equal(pipe$stats$unmapped, 0)
  expect_identical(pipe$multiqc_status, "unavailable")
  expect_true(file.exists(paste0(pipe$bam_files, ".bai")))
  expect_false(any(grepl("PEAR|/pear", pipe$commands)))
  cli_expect_extracted(pipe, f)
})

test_that("mixed tables preserve paired and R1-only samples including indel barcodes", {
  cli_tools_available(paired = TRUE)
  f <- cli_read_fixture(); on.exit(unlink(f$folder, recursive = TRUE))
  lookup <- .barbac_tool
  local_mocked_bindings(.barbac_tool = function(tool, ...) if (tool == "multiqc") "" else lookup(tool, ...))
  samples <- data.frame(sample = c("paired", "single_na", "single_empty"),
    R1 = f$r1, R2 = c(f$r2, NA, ""))
  pipe <- run_cli_pipeline(samples, f$reference, file.path(f$folder, "mixed"), verbose = FALSE)
  expect_identical(names(pipe$bam_files), samples$sample)
  expect_identical(pipe$samples$mode, c("paired", "single", "single"))
  expect_equal(pipe$stats$mapped, rep(144, 3))
  expect_equal(pipe$stats$unmapped, rep(0, 3))
  expect_true(file.exists(file.path(pipe$merged_dir, "paired_ANC.assembled.fastq")))
  expect_true(file.exists(file.path(pipe$bam_dir, "paired_ANC.assembled_sorted.bam")))
  cli_expect_extracted(pipe, f)
})

test_that("invalid inputs and existing outputs fail before starting external work", {
  f <- cli_read_fixture(); on.exit(unlink(f$folder, recursive = TRUE))
  samples <- data.frame(sample = "one", R1 = f$r1)
  out <- file.path(f$folder, "results"); dir.create(out)
  writeLines("preserve me", file.path(out, "old.txt"))
  expect_error(run_cli_pipeline(samples, f$reference, out), "new or empty")
  expect_identical(readLines(file.path(out, "old.txt")), "preserve me")
  expect_error(run_cli_pipeline(rbind(samples, samples), f$reference), "unique")
  samples$sample <- "../outside"
  expect_error(run_cli_pipeline(samples, f$reference), "Sample names")
  samples$sample <- "one"; samples$R1 <- "missing.fastq"
  expect_error(run_cli_pipeline(samples, f$reference), "existing FASTQ")
  expect_error(run_cli_pipeline(data.frame(sample = character(), R1 = character()), f$reference), "at least one row")
})

test_that("a failed command stops the pipeline and records its diagnostics", {
  skip_on_os("windows")
  f <- cli_read_fixture(); on.exit(unlink(f$folder, recursive = TRUE))
  executable <- file.path(f$folder, "failing-tool")
  writeLines(c("#!/bin/sh", "echo 'deliberate failure' >&2", "exit 17"), executable)
  Sys.chmod(executable, "0755")
  local_mocked_bindings(.barbac_require_tools = function(tools, ...) setNames(rep(executable, length(tools)), tools),
                        .barbac_tool = function(...) "")
  out <- file.path(f$folder, "failed")
  expect_error(run_cli_pipeline(data.frame(sample = "one", R1 = f$r1), f$reference, out, verbose = FALSE),
               "fastqc failed with exit code 17")
  expect_match(paste(readLines(file.path(out, "pipeline.log")), collapse = "\n"), "deliberate failure")
  expect_length(list.files(file.path(out, "merged", "bam")), 0)
})
