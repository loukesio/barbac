test_that("barbac_ts_area returns a ggplot for long-format input", {
  df <- data.frame(
    barcode = rep(c("A", "B", "C"), each = 3),
    time    = rep(0:2, times = 3),
    counts  = c(100, 80, 60,
                 20, 40, 60,
                  0,  0, 30)
  )
  p <- barbac_ts_area(df, min_total_count = 0,
                      palette = c("red", "blue", "green"))
  expect_s3_class(p, "ggplot")
  # frequencies should sum to ~1 at each timepoint we haven't zeroed out
  totals <- tapply(p$data$.freq, p$data$time, sum)
  expect_true(all(abs(totals - 1) < 1e-5 | totals < 1e-5))
})

test_that("barbac_ts_area auto-detects Bartender wide format", {
  df <- data.frame(
    Cluster.ID    = c("bc1", "bc2"),
    Center        = c("ACGT", "ACGA"),
    Cluster.Score = c(0.9, 0.8),
    time_point_1  = c(100, 20),
    time_point_2  = c( 80, 40),
    time_point_3  = c( 60, 60)
  )
  p <- barbac_ts_area(df, min_total_count = 0,
                      palette = c("black", "grey"))
  expect_s3_class(p, "ggplot")
  # With time_zero_shift = FALSE (default), the wide time_point_1
  # column keeps its native value of 1.
  expect_true(1 %in% p$data$time)

  # With time_zero_shift = TRUE, time == 1 gets remapped to 0.
  p_shift <- barbac_ts_area(df, min_total_count = 0,
                            palette = c("black", "grey"),
                            time_zero_shift = TRUE)
  expect_true(0 %in% p_shift$data$time)
})

test_that("barbac_ts_area errors informatively on unknown schema", {
  bad <- data.frame(foo = 1:3, bar = 4:6)
  expect_error(barbac_ts_area(bad),
               regexp = "Unrecognised input layout")
})

test_that("include_late = FALSE drops lineages absent at t = 0", {
  # A is observed at t = 0, 1, 2. B appears only at t = 2.
  df <- data.frame(
    barcode = c("A", "A", "A", "B"),
    time    = c(0,   1,   2,   2),
    counts  = c(100, 80,  60,  90)
  )
  p <- barbac_ts_area(df, min_total_count = 0,
                      include_late = FALSE,
                      time_zero_shift = FALSE,
                      palette = c("red", "blue"))
  expect_setequal(unique(as.character(p$data$barcode)), "A")
})

test_that("fill_missing = 'zero' fills late-arriving cells with 0", {
  df <- data.frame(
    barcode = c("A", "B"),
    time    = c(2,  2),          # both only observed at t = 2
    counts  = c(50, 50)
  )
  df <- rbind(df, data.frame(barcode = "A", time = 0, counts = 100))
  p <- barbac_ts_area(df, min_total_count = 0,
                      include_late = TRUE,
                      fill_missing = "zero",
                      time_zero_shift = FALSE,
                      palette = c("red", "blue"))
  # B at time 0 should be filled with 0 exactly, not epsilon
  b_at_zero <- p$data$.freq[p$data$barcode == "B" & p$data$time == 0]
  expect_equal(b_at_zero, 0)
})
