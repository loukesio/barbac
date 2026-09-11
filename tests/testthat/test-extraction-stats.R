test_that("extraction statistics bin lengths by requested bounds, independent of row order", {
  data <- data.frame(barcode = c("AAAA", "AAA", "AAAAAA", "AAAAA", "AAAAAAA"),
                     counts = c(2, 1, 4, 3, 5), barcode_length = c(4, 3, 6, 5, 7))
  output <- capture.output(plot <- barbac_xtr.stats(data, c(4, 6), verbose = TRUE))
  expect_s3_class(plot, "patchwork")
  expect_match(paste(output, collapse = "\n"), "4 - 6 bp\\s+3\\s+9")
  expect_match(paste(output, collapse = "\n"), "< 4 bp\\s+1\\s+1")
  expect_match(paste(output, collapse = "\n"), "> 6 bp\\s+1\\s+5")
  reordered <- capture.output(reordered_plot <- barbac_xtr.stats(data[c(5, 3, 1, 4, 2), ], c(4, 6), verbose = TRUE))
  expect_identical(output, reordered)
})

test_that("report details distinguish sequence entropy from abundance and retain counts", {
  data <- data.frame(barcode = c("AAAA", "ACGT", "ACAC"),
    counts = c(100, 1, 9), barcode_length = c(4, 4, 4))
  result <- barbac_xtr.stats(data, c(3, 5), verbose = FALSE,
                            panel_labels = TRUE, return_details = TRUE)
  expect_s3_class(result$plot, "patchwork")
  expect_equal(unname(result$sequence_entropy$entropy), c(0, 2, 1))
  expect_equal(result$length_summary$barcodes, 3)
  expect_equal(result$length_summary$reads, 110)
  expect_named(result$plots, c("length", "abundance", "entropy"))
})
