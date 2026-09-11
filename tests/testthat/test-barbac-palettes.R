test_that("all 32 LTC palettes are available without another palette package", {
  palettes <- barbac_palettes()
  expect_length(palettes, 32)
  expect_false(anyDuplicated(names(palettes)) > 0)
  expect_identical(palettes$alger,
                   c("#1A5B5B", "#ACC8BE", "#F4AB5C", "#D1422F"))
  for (name in names(palettes)) {
    colours <- barbac_palette(name, 100)
    expect_length(colours, 100)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colours)))
    expect_equal(grDevices::col2rgb(colours[c(1, 100)]),
                 grDevices::col2rgb(palettes[[name]][c(1, length(palettes[[name]]))]))
  }
  palettes$alger[1] <- "red"
  expect_identical(barbac_palette("alger")[1], "#1A5B5B")
})

test_that("palette names accept aliases and interpolation covers both ends", {
  for (alias in c("Casa Natal", "CASA_NATAL", "casa-natal", "casanatal")) {
    expect_identical(barbac_palette(alias, 9), barbac_palette("casa_natal", 9))
  }
  expect_identical(barbac_palette("alger", 2), c("#1A5B5B", "#D1422F"))
  expect_identical(barbac_palette("alger", 1), "#1A5B5B")
  expect_identical(barbac_palette("alger", 0), character())
  expect_error(barbac_palette("algre"), "Unknown LTC palette")
  for (bad_name in list(NA_character_, character(), c("alger", "dora"), 1)) {
    expect_error(barbac_palette(bad_name), "single LTC palette name")
  }
  for (n in list(-1, 1.5, Inf, NA_real_, numeric(), c(1, 2), "2")) {
    expect_error(barbac_palette("alger", n), "non-negative whole number")
  }
})

test_that("named palettes preserve frequencies and individual barcode groups", {
  df <- data.frame(barcode = rep(c("b", "a", "c"), each = 3),
                   time = rep(0:2, 3), counts = c(9, 6, 3, 1, 3, 5, 0, 1, 2))
  reference <- barbac_ts_area(df, min_total_count = 0, fill_missing = "zero",
                              palette = c("red", "green", "blue"), theme = NULL)
  for (name in names(barbac_palettes())) {
    p <- barbac_ts_area(df, min_total_count = 0, fill_missing = "zero",
                        palette = name, theme = NULL)
    expect_identical(p$data, reference$data)
    expect_identical(unname(p$scales$get_scales("fill")$palette(3)),
                     barbac_palette(name, 3))
  }
  built <- ggplot2::ggplot_build(p)$data[[1]]
  expect_length(unique(built$group), 3)
  expect_equal(max(built$ymax), 1)
})

test_that("existing default, custom vectors and single colours still work", {
  df <- data.frame(barcode = c("a", "b", "c"), time = 0, counts = 1)
  colours <- function(palette = NULL) {
    p <- barbac_ts_area(df, min_total_count = 0, palette = palette, theme = NULL)
    unname(p$scales$get_scales("fill")$palette(3))
  }
  base <- PNWColors::pnw_palette("Sailboat", n = 3, type = "continuous")
  expect_identical(colours(), grDevices::colorRampPalette(base)(3))
  expect_identical(colours(c("red", "blue", "green", "black")),
                   c("red", "blue", "green"))
  expect_identical(colours(c("black", "white")), c("#000000", "#7F7F7F", "#FFFFFF"))
  expect_identical(colours("red"), rep("#FF0000", 3))
  expect_identical(colours("#123456"), rep("#123456", 3))
  expect_error(colours("algre"), "Unknown palette name or invalid colour")
  expect_error(colours(c("red", "alger")), "Unknown palette name or invalid colour")
  expect_error(colours(character()), "non-empty colour vector")
  expect_error(colours(NA_character_), "non-empty colour vector")
})
