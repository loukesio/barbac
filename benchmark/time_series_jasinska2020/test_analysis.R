#!/usr/bin/env Rscript
source('benchmark/time_series_jasinska2020/analysis_helpers.R')
source('benchmark/time_series_jasinska2020/frequency_audit.R')
testthat::test_that('diversity uses positive counts and explicit denominators', {
  d <- diversity(c(25,25,25,25,0),100)
  testthat::expect_equal(unname(d),rep(4,5))
  testthat::expect_equal(unname(diversity(c(100,0),100)),rep(1,5))
  d <- diversity(c(10,10),100)
  testthat::expect_equal(unname(d['assigned_shannon_effective']),2)
  testthat::expect_equal(unname(d['inverse_dominance']),10)
  testthat::expect_error(diversity(c(100,1),100))
  testthat::expect_error(diversity(c(-1,3),100))
})
testthat::test_that('publication reference retains short barcode identities', {
  a <- read_author_tables('benchmark/time_series_jasinska2020')
  testthat::expect_equal(nrow(a$top),300)
  testthat::expect_true('AGCAGAAAGGTGC' %in% a$top$barcode)
  testthat::expect_true('TCAATTGCCACCAA' %in% a$top$barcode)
})
testthat::test_that('frequency audit distinguishes denominators and retains absent calls', {
  top <- data.frame(population='P',barcode=c('AAAA','CCCC'),
    published_final=c(.6,.3),barbac_final=c(.48,.24))
  samples <- data.frame(population='P',passage=30,input_reads_barbac=100,
    extracted_reads=80,extracted_fraction=.8)
  r <- compare_frequency_denominators(top,samples)
  testthat::expect_equal(r$top$barbac_final_extracted,c(.6,.3))
  testthat::expect_equal(r$summary$Mean_absolute_difference_input_pp,9)
  testthat::expect_equal(r$summary$Mean_absolute_difference_extracted_pp,0,tolerance=1e-12)
  testthat::expect_true(r$audit$Incompatible_with_all_input)
  top$barbac_final[2] <- 0
  r <- compare_frequency_denominators(top,samples)
  testthat::expect_equal(nrow(r$top),2)
  testthat::expect_equal(r$top$barbac_final_extracted[2],0)
  top$published_final <- c(.4,.40004)
  testthat::expect_false(compare_frequency_denominators(top,samples)$audit$Incompatible_with_all_input)
  top$published_final[2] <- .40006
  testthat::expect_true(compare_frequency_denominators(top,samples)$audit$Incompatible_with_all_input)
})
# Provenance must fail closed if JSON drops filename keys or an input changes.
source('benchmark/time_series_jasinska2020/provenance_helpers.R')
testthat::test_that('checksum maps survive JSON and reject missing keys or changed files', {
  dir <- tempfile('barbac-provenance-');dir.create(dir)
  file <- file.path(dir,'counts.csv');writeLines('barcode,counts\nACGT,7',file)
  hashes <- list('counts.csv'=digest::digest(file=file,algo='sha256'))
  json <- file.path(dir,'hashes.json')
  jsonlite::write_json(hashes,json,auto_unbox=TRUE)
  decoded <- jsonlite::read_json(json,simplifyVector=TRUE)
  testthat::expect_true(assert_hash_map(decoded,dir))
  testthat::expect_error(assert_hash_map(unname(unlist(decoded)),dir))
  testthat::expect_error(assert_hash_map(list(),dir))
  writeLines('barcode,counts\nACGT,8',file)
  testthat::expect_error(assert_hash_map(decoded,dir))
})
