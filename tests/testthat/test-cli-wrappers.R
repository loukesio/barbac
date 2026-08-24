test_that("run_multiqc quotes paths and reports command failures", {
  input_dir <- file.path(tempdir(), "fastqc reports")
  output_dir <- file.path(tempdir(), "multiqc reports")
  dir.create(input_dir, showWarnings = FALSE)
  seen <- NULL

  local_mocked_bindings(
    .barbac_which = function(x) setNames("/mock/multiqc", x),
    .barbac_system = function(command, ...) {
      seen <<- command
      0L
    },
    .package = "barbac"
  )
  run_multiqc(input_dir, output_dir)

  expect_identical(
    seen,
    sprintf("multiqc %s -o %s", shQuote(input_dir), shQuote(output_dir))
  )

  local_mocked_bindings(
    .barbac_which = function(x) setNames("/mock/multiqc", x),
    .barbac_system = function(...) 2L,
    .package = "barbac"
  )
  expect_error(run_multiqc(input_dir, output_dir), "exit code 2")
})

test_that("run_fastqc quotes FASTQ and output paths", {
  base_dir <- file.path(tempdir(), "fastqc input directory")
  dir.create(base_dir, showWarnings = FALSE)
  fastq <- file.path(base_dir, "sample one.fastq")
  samples <- file.path(base_dir, "samples.csv")
  writeLines(c("@read", "ACGT", "+", "IIII"), fastq)
  write.csv(data.frame(R1 = fastq), samples, row.names = FALSE)
  seen <- NULL

  local_mocked_bindings(
    .barbac_which = function(x) setNames("/mock/fastqc", x),
    .barbac_system = function(command, ...) {
      seen <<- command
      0L
    },
    .package = "barbac"
  )
  run_fastqc(samples)

  expected_output <- file.path(base_dir, "results", "fastQC")
  expect_identical(
    seen,
    sprintf("fastqc %s -o %s", shQuote(fastq), shQuote(expected_output))
  )
})
