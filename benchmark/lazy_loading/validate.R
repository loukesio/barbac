args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
.libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages(library(barbac))
stopifnot(!'GenomicAlignments' %in% loadedNamespaces())
stopifnot(!'GenomicRanges' %in% loadedNamespaces())
x <- as.data.frame(testthat::test_dir('tests/testthat', package = 'barbac', reporter = 'summary', stop_on_failure = TRUE))
receipt <- list(cases = nrow(x), assertions = sum(x$passed), failures = sum(x$failed),
  errors = sum(x$error), skipped = sum(x$skipped), warnings = sum(x$warning),
  native_build = barbac:::barbac_build_id(), library = normalizePath(args[[1]]))
stopifnot(receipt$failures == 0, receipt$errors == 0, receipt$skipped == 0)
jsonlite::write_json(receipt, 'benchmark/lazy_loading/package_tests.json', auto_unbox = TRUE, pretty = TRUE)
