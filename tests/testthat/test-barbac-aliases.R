test_that("barbac_cluster() is identical to super_cluster2()", {
  input <- data.frame(
    barcode = c("AAAAAAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAAAAAAT",
                "AAAAAAAAAAAAAAAAAATT",
                "TTTTTTTTTTTTTTTTTTTT",
                "TTTTTTTTTTTTTTTTTTTA",
                "GGGGGGGGGGGGGGGGGGG"),   # one deletion from a G20 parent
    counts  = c(1000, 10, 3, 800, 8, 5)
  )
  for (method in c("lv", "hamming")) {
    # Hamming warns about the mixed-length input; the warning is expected and
    # must be identical for both names, so it is silenced here.
    original <- suppressWarnings(super_cluster2(input, distance = 3, method = method, verbose = FALSE))
    alias <- suppressWarnings(barbac_cluster(input, distance = 3, method = method, verbose = FALSE))
    expect_identical(alias, original)
  }
  original <- super_cluster2(input, distance = 2, merge_ratio = 5,
                             tie_break = "support", indel_model = "poisson",
                             verbose = FALSE)
  alias <- barbac_cluster(input, distance = 2, merge_ratio = 5,
                          tie_break = "support", indel_model = "poisson",
                          verbose = FALSE)
  expect_identical(alias, original)
})

test_that("aliases forward to the original functions", {
  expect_identical(formals(barbac_cluster), formals(function(...) NULL))
  expect_identical(formals(barbac_extract), formals(function(...) NULL))
  expect_identical(formals(barbac_extract_stats), formals(function(...) NULL))
  expect_identical(formals(barbac_map), formals(function(...) NULL))
  expect_identical(body(barbac_cluster), quote({ super_cluster2(...) }))
  expect_identical(body(barbac_extract), quote({ barbac_xtr(...) }))
  expect_identical(body(barbac_extract_stats), quote({ barbac_xtr.stats(...) }))
  expect_identical(body(barbac_map), quote({ run_cli_pipeline(...) }))
})

test_that("barbac_extract_stats() matches barbac_xtr.stats() on a count table", {
  counts <- data.frame(
    barcode = c("ACGTACGTAC", "ACGTACGTAG", "ACGTACGTACG", "ACGTACGT"),
    counts = c(100L, 20L, 5L, 2L),
    barcode_length = c(10L, 10L, 11L, 8L)
  )
  original <- barbac_xtr.stats(counts, barcode_length = c(9, 11), verbose = FALSE)
  alias <- barbac_extract_stats(counts, barcode_length = c(9, 11), verbose = FALSE)
  expect_identical(alias$length_summary, original$length_summary)
  expect_identical(alias$sequence_entropy, original$sequence_entropy)
})
