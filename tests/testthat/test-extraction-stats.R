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
