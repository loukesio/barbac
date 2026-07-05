test_that("super_cluster2 returns the expected result shape", {
  input <- data.frame(
    barcode = c("AAAAAAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAAAAAAT",   # one sub
                "TTTTTTTTTTTTTTTTTTTT",
                "TTTTTTTTTTTTTTTTTTTA"),  # one sub
    counts  = c(1000, 10, 800, 8)
  )
  res <- super_cluster2(input, distance = 3, verbose = FALSE)

  expect_s3_class(res, "tbl_df")
  expect_true(all(c("cluster_id", "central_barcode",
                    "all_barcodes", "all_counts", "sum_counts") %in% names(res)))
  expect_type(res$all_barcodes, "list")
  expect_type(res$all_counts,   "list")
})

test_that("super_cluster2 with distance = 0 keeps every unique input", {
  input <- data.frame(
    barcode = c("AAAA", "AAAT", "AATA", "TTTT"),
    counts  = c(100, 50, 25, 10)
  )
  res <- super_cluster2(input, distance = 0, verbose = FALSE)

  expect_equal(nrow(res), 4L)
  expect_setequal(res$central_barcode, input$barcode)
})

test_that("super_cluster2 with distance >= 1 collapses single-substitution neighbours", {
  input <- data.frame(
    barcode = c("AAAAAAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAAAAAAT"),  # exactly 1 sub away
    counts  = c(1000, 10)                 # child abundance is dominated
  )
  res <- super_cluster2(input, distance = 3, verbose = FALSE)

  expect_equal(nrow(res), 1L)
  expect_equal(res$central_barcode, "AAAAAAAAAAAAAAAAAAAA")
  expect_equal(res$sum_counts, 1010)
})

test_that("super_cluster2 preserves total counts", {
  set.seed(42)
  input <- data.frame(
    barcode = replicate(20, paste0(sample(c("A","C","G","T"), 20, replace = TRUE),
                                   collapse = "")),
    counts  = sample(1:1000, 20)
  )
  res <- super_cluster2(input, distance = 3, verbose = FALSE)
  expect_equal(sum(res$sum_counts), sum(input$counts))
})

test_that("super_cluster2 rejects unsupported distance methods", {
  input <- data.frame(barcode = c("AAAA", "AAAT"), counts = c(10, 5))
  expect_error(super_cluster2(input, method = "jw",     verbose = FALSE))
  expect_error(super_cluster2(input, method = "cosine", verbose = FALSE))
  expect_error(super_cluster2(input, method = "qgram",  verbose = FALSE))
})

test_that("super_cluster2 collapses exact-duplicate barcodes into one cluster", {
  input <- data.frame(
    barcode = c("AAAAAAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAAAAAAA",   # identical duplicate
                "TTTTTTTTTTTTTTTTTTTT"),
    counts  = c(100, 50, 30)
  )
  res <- super_cluster2(input, distance = 0, verbose = FALSE)

  expect_equal(nrow(res), 2L)
  a <- res[res$central_barcode == "AAAAAAAAAAAAAAAAAAAA", ]
  expect_equal(a$sum_counts, 150L)              # 100 + 50 pooled, not double-counted
  expect_equal(sum(res$sum_counts), 180L)
})

test_that("super_cluster2 warns on Hamming-incompatible barcodes", {
  input <- data.frame(
    barcode = c("AAAAAAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAAAAAAN"),   # contains N, not 2-bit packable
    counts  = c(100, 5)
  )
  expect_warning(super_cluster2(input, method = "hamming", verbose = FALSE),
                 "Hamming mode cannot compare")
})
