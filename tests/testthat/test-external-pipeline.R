test_that("real single-end FASTQ reaches mapping and BAM summary", {
  skip_if(tolower(Sys.getenv("BARBAC_RUN_EXTERNAL_TESTS")) != "true",
          "opt-in external-tool integration test")
  needed <- c("fastqc", "minimap2", "samtools")
  skip_if(any(Sys.which(needed) == ""), "required external tools unavailable")

  work <- tempfile("barbac external ")
  dir.create(work)
  reference <- file.path(work, "reference sequence.fa")
  fastq <- file.path(work, "single read.fastq")
  writeLines(c(">ref", "ACGTACGTACGTACGTACGT"), reference)
  writeLines(c("@read1", "ACGTACGTACGTACGTACGT", "+", "IIIIIIIIIIIIIIIIIIII"), fastq)

  result <- run_cli_pipeline(
    data.frame(sample = "sample one", R1 = fastq),
    reference = reference,
    output_dir = file.path(work, "pipeline output"),
    verbose = FALSE
  )

  expect_true(file.exists(result$summary_file))
  expect_equal(nrow(result$stats), 1L)
  expect_equal(result$stats$mapped + result$stats$unmapped, 1L)
})
