test_that("single-end R1 is routed directly to minimap2", {
  merged_dir <- tempfile("barbac-merged-")
  dir.create(merged_dir)
  r1 <- tempfile(fileext = ".fastq")
  writeLines(c("@read", "ACGT", "+", "IIII"), r1)

  samples <- data.frame(sample = "single_sample", R1 = r1)
  inputs <- barbac:::.pipeline_mapping_inputs(samples, merged_dir)

  expect_identical(inputs$files, r1)
  expect_identical(inputs$names, "single_sample")
  expect_length(inputs$merged_files, 0L)
})

test_that("paired-end samples still use PEAR assembled output", {
  merged_dir <- tempfile("barbac-merged-")
  dir.create(merged_dir)
  assembled <- file.path(merged_dir, "paired_ANC.assembled.fastq")
  writeLines(c("@read", "ACGT", "+", "IIII"), assembled)

  samples <- data.frame(
    sample = "paired",
    R1 = "paired_R1.fastq",
    R2 = "paired_R2.fastq"
  )
  inputs <- barbac:::.pipeline_mapping_inputs(samples, merged_dir)

  expect_identical(inputs$files, assembled)
  expect_identical(inputs$names, "paired_ANC.assembled")
  expect_length(inputs$single_files, 0L)
})

test_that("mixed input routes paired assemblies and single-end R1", {
  merged_dir <- tempfile("barbac-merged-")
  dir.create(merged_dir)
  assembled <- file.path(merged_dir, "paired_ANC.assembled.fastq")
  single_r1 <- tempfile(fileext = ".fastq")
  writeLines(c("@read", "ACGT", "+", "IIII"), assembled)
  writeLines(c("@read", "TGCA", "+", "IIII"), single_r1)

  samples <- data.frame(
    sample = c("paired", "single"),
    R1 = c("paired_R1.fastq", single_r1),
    R2 = c("paired_R2.fastq", NA_character_)
  )
  inputs <- barbac:::.pipeline_mapping_inputs(samples, merged_dir)

  expect_identical(inputs$files, c(assembled, single_r1))
  expect_identical(inputs$names, c("paired_ANC.assembled", "single"))
})
