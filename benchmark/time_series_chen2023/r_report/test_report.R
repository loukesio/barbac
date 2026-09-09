# Run after loading the project package; no sequencing files needed.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
here <- dirname(normalizePath(script))
devtools::load_all(normalizePath(file.path(here, "..", "..", "..")), quiet = TRUE)
source(file.path(here, "report_helpers.R"))
testthat::test_that("zero observations do not become observed clusters", {
  s <- positive_cluster_stats(c(0, 0, 1, 9))
  testthat::expect_equal(s$n_clusters, 2)
  testthat::expect_equal(s$total_reads, 10)
  testthat::expect_equal(s$n_singletons, 1)
  testthat::expect_error(positive_cluster_stats(c(1, -1)))
})
testthat::test_that("Other conserves the denominator and zero timepoints remain zero", {
  d <- data.frame(Barcode = rep(c("A", "B", "C"), 2), generation = rep(c(8, 40), each = 3),
                  counts = c(10, 0, 0, 0, 2, 8))
  x <- composition_input(d, c("A", "B"), c(A = "L01", B = "L02"))
  testthat::expect_equal(as.numeric(tapply(x$counts, x$time, sum)), c(10, 10))
  testthat::expect_equal(x$counts[x$barcode == "Other" & x$time == 40], 8)
  p <- barbac::barbac_ts_area(x, min_total_count = 0, fill_missing = "zero",
    time_zero_shift = FALSE, palette = c("#2C4A63", "#B07D22", "#D4D4D4"))
  testthat::expect_equal(p$data$.freq[p$data$barcode == "L02" & p$data$time == 8], 0)
  testthat::expect_equal(p$data$.freq[p$data$barcode == "L02" & p$data$time == 40], .2)
  testthat::expect_setequal(p$data$time, c(8, 40))
  testthat::expect_error(composition_input(rbind(d, d[1, ]), c("A", "B"), c(A = "L01", B = "L02")))
})
