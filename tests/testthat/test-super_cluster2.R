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

test_that("super_cluster2 output is invariant to input row order", {
  set.seed(3)
  ALPH  <- c("A", "C", "G", "T")
  truth <- unique(replicate(120, paste0(sample(ALPH, 20, TRUE), collapse = "")))
  mk <- function(s) {                       # one random substitution
    cs <- strsplit(s, "")[[1]]; j <- sample(20, 1)
    cs[j] <- sample(setdiff(ALPH, cs[j]), 1); paste0(cs, collapse = "")
  }
  bc  <- unlist(lapply(truth, function(t) c(t, replicate(2, mk(t)))))
  tab <- aggregate(counts ~ barcode,
                   data.frame(barcode = bc,
                              counts  = sample(1:5, length(bc), TRUE)),
                   sum)
  canon <- function(d) {
    r <- super_cluster2(d, distance = 3, verbose = FALSE)
    sort(vapply(seq_len(nrow(r)), function(i)
      paste(r$central_barcode[i], r$sum_counts[i],
            paste(sort(r$all_barcodes[[i]]), collapse = "|"), sep = "~"),
      character(1)))
  }
  base <- canon(tab)
  expect_identical(base, canon(tab[sample(nrow(tab)), ]))   # rows shuffled
  expect_identical(base, canon(tab[order(tab$counts), ]))   # count-ascending
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

test_that("indexed LV clustering matches the full scan on fixed-anchor designs", {
  # Structured barcode designs repeat a constant anchor in every sequence, so
  # seeds drawn from the anchor are shared by the whole table and carry no
  # information. The seed index skips such posting lists to stay fast; this
  # checks that skipping them costs no recall by comparing against the
  # exhaustive scan (use_kmer_filter = FALSE), which uses no index at all.
  set.seed(11)
  ALPH <- c("A", "C", "G", "T")
  template <- c(rep("N", 8), strsplit("ATGC", "")[[1]],
                rep("N", 8), strsplit("ATCGTTAA", "")[[1]])
  var_pos <- which(template == "N")

  draw <- function() {
    s <- template
    s[var_pos] <- sample(ALPH, length(var_pos), replace = TRUE)
    paste0(s, collapse = "")
  }
  truth <- unique(replicate(60, draw()))

  mutate1 <- function(s) {                  # one substitution in the variable part
    cs <- strsplit(s, "")[[1]]
    j <- sample(var_pos, 1)
    cs[j] <- sample(setdiff(ALPH, cs[j]), 1)
    paste0(cs, collapse = "")
  }
  bc  <- unlist(lapply(truth, function(t) c(t, t, mutate1(t))))
  tab <- aggregate(counts ~ barcode,
                   data.frame(barcode = bc,
                              counts  = sample(1:40, length(bc), TRUE)),
                   sum)

  canon <- function(filter_on) {
    r <- super_cluster2(tab, distance = 3, method = "lv",
                        use_kmer_filter = filter_on, verbose = FALSE)
    sort(vapply(seq_len(nrow(r)), function(i)
      paste(r$central_barcode[i], r$sum_counts[i],
            paste(sort(r$all_barcodes[[i]]), collapse = "|"), sep = "~"),
      character(1)))
  }

  expect_identical(canon(TRUE), canon(FALSE))
})

test_that("tie_break = 'hash' stays deterministic and order-invariant", {
  set.seed(7)
  ALPH <- c("A", "C", "G", "T")
  truth <- replicate(40, paste0(sample(ALPH, 20, replace = TRUE), collapse = ""))
  mk <- function(s) {
    cs <- strsplit(s, "")[[1]]; j <- sample(20, 1)
    cs[j] <- sample(setdiff(ALPH, cs[j]), 1); paste0(cs, collapse = "")
  }
  bc  <- unlist(lapply(truth, function(t) c(t, mk(t))))
  tab <- aggregate(counts ~ barcode,
                   data.frame(barcode = bc,
                              counts  = sample(1:6, length(bc), TRUE)),
                   sum)

  canon <- function(d, ...) {
    r <- super_cluster2(d, distance = 3, verbose = FALSE, ...)
    sort(paste(r$central_barcode, r$sum_counts, sep = "~"))
  }

  # A seed is reproducible, and still independent of input row order.
  a <- canon(tab,                      tie_break = "hash", tie_seed = 3L)
  expect_identical(a, canon(tab,       tie_break = "hash", tie_seed = 3L))
  expect_identical(a, canon(tab[sample(nrow(tab)), ],
                                       tie_break = "hash", tie_seed = 3L))

  # The default is unchanged by the option's existence.
  expect_identical(canon(tab), canon(tab, tie_break = "sequence"))

  # Seeds are genuinely different orderings, not the same one relabelled.
  keys <- barbac:::barbac_seq_order_key(tab$barcode, 1L)
  expect_false(identical(keys, barbac:::barbac_seq_order_key(tab$barcode, 2L)))
  expect_identical(keys, barbac:::barbac_seq_order_key(tab$barcode, 1L))
})

test_that("Hamming mode absorbs trace indel reads instead of splitting them off", {
  # The Hamming partition index is keyed by sequence length, so a read carrying
  # an indel is never offered its parent as a candidate and would found a
  # cluster of its own. With a trace of such reads they are rescued by edit
  # distance and land in the right cluster.
  set.seed(3)
  ALPH <- c("A", "C", "G", "T")
  truth <- unique(replicate(120, paste0(sample(ALPH, 20, replace = TRUE),
                                        collapse = "")))
  parent <- truth[1]
  deletion <- substr(parent, 1, 19)          # one base short: an indel read

  input <- data.frame(
    barcode = c(truth, deletion),
    counts  = c(rep(200L, length(truth)), 1L),
    stringsAsFactors = FALSE
  )

  res <- super_cluster2(input, distance = 3, method = "hamming", verbose = FALSE)

  # The shortened read must not survive as its own centroid...
  expect_false(deletion %in% res$central_barcode)
  # ...it belongs to the barcode it was derived from.
  owner <- res$all_barcodes[[which(res$central_barcode == parent)]]
  expect_true(deletion %in% owner)
})

test_that("Hamming mode warns when the data is substantially length-variable", {
  # Past a trace, rescuing by edit distance is the wrong answer: the data wants
  # Levenshtein, and the user should be told so rather than handed a slow run.
  set.seed(4)
  ALPH <- c("A", "C", "G", "T")
  truth <- unique(replicate(40, paste0(sample(ALPH, 20, replace = TRUE),
                                       collapse = "")))
  shortened <- substr(truth, 1, 19)          # half the table is length 19
  input <- data.frame(
    barcode = c(truth, shortened),
    counts  = c(rep(50L, length(truth)), rep(3L, length(shortened))),
    stringsAsFactors = FALSE
  )
  expect_warning(super_cluster2(input, distance = 3, method = "hamming",
                                verbose = FALSE),
                 "differ from the modal barcode length")
})
