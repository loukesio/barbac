canonical_members <- function(x) {
  sort(vapply(seq_len(nrow(x)), function(i)
    paste(x$central_barcode[i], x$sum_counts[i],
          paste(sort(x$all_barcodes[[i]]), collapse = '|')), character(1)))
}

test_that('Poisson guard responds to error rate and repeated gap opportunities', {
  # Synthetic sequences and counts fixed independently of the reference cases.
  parent <- 'GTAAAAAAAAAACTG'
  child <- 'GTAAAAAAAAACTG'
  input <- data.frame(barcode = c(parent, child), counts = c(1000L, 65L))
  cluster <- function(model, rate = 0.005, tab = input) super_cluster2(
    tab, method = 'lv', distance = 1, indel_model = model,
    error_rate = rate, verbose = FALSE)
  expect_equal(nrow(cluster('none')), 2)
  expect_equal(cluster('poisson')$central_barcode, parent)
  expect_equal(cluster('poisson')$sum_counts, 1065L)
  expect_equal(nrow(cluster('poisson', 0.0005)), 2)
  # A nearby well-supported true length variant must remain separate.
  input$counts[2] <- 300L
  expect_equal(nrow(cluster('poisson', tab = input)), 2)
  # Ordinary substitutions receive no relaxation from this option.
  input$barcode[2] <- 'GTAAAAACAAAACTG'
  input$counts[2] <- 65L
  expect_identical(cluster('none', tab = input), cluster('poisson', tab = input))
})

test_that('single-base insertion and boundary-run deletion models are directional', {
  # Specific insertions have a quarter of the deletion event rate.
  tab <- data.frame(barcode = c('GTAAAAAAAAAACTG', 'GTAAAAAAAAAAACTG'),
                    counts = c(1000L, 65L))
  r <- super_cluster2(tab, method = 'lv', distance = 1,
                      indel_model = 'poisson', verbose = FALSE)
  expect_equal(nrow(r), 2)
  for (parent in c('AAAAAAAAAAGT', 'GTAAAAAAAAAA', 'AAAAAAAAAA')) {
    child <- sub('A', '', parent)
    tab <- data.frame(barcode = c(parent, child), counts = c(1000L, 65L))
    r <- super_cluster2(tab, distance = 1, indel_model = 'poisson', verbose = FALSE)
    expect_equal(r$central_barcode, parent)
    expect_equal(r$sum_counts, 1065L)
  }
})

test_that('Poisson option is explicit and propagates through files and directories', {
  tab <- data.frame(barcode = c('GTAAAAAAAAAACTG', 'GTAAAAAAAAACTG'), counts = c(1000L,65L))
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  path <- file.path(root, 'input.csv'); write.csv(tab, path, row.names = FALSE)
  direct <- super_cluster2(tab, indel_model = 'poisson', verbose = FALSE)
  from_file <- super_cluster2(path, indel_model = 'poisson', verbose = FALSE)
  from_dir <- super_cluster2(root, indel_model = 'poisson', verbose = FALSE)[[1]]
  expect_identical(canonical_members(direct), canonical_members(from_file))
  expect_identical(canonical_members(direct), canonical_members(from_dir))
  expect_error(super_cluster2(tab, indel_model = 'bad'), 'arg')
  expect_error(super_cluster2(tab, method = 'hamming', indel_model = 'poisson'), 'requires')
  expect_error(super_cluster2(tab, use_cpp = FALSE, indel_model = 'poisson'), 'requires')
})

test_that('abundance-pruned LV agrees with full scans with either indel model', {
  set.seed(7249)
  roots <- replicate(35, paste(sample(c('A','C','G','T'), 20, TRUE), collapse = ''))
  roots <- c(roots, 'GTAAAAAAAAAACTG', 'GTAAAAAAAAACTG')
  variants <- unlist(lapply(roots, function(s) c(s,
    paste0(substr(s, 1, 4), substr(s, 6, nchar(s))),
    paste0(substr(s, 1, 6), 'A', substr(s, 7, nchar(s))),
    paste0('G', substr(s, 2, nchar(s))))))
  tab <- data.frame(barcode = variants, counts = sample(c(1L,2L,10L,65L,1000L,1000000L), length(variants), TRUE))
  for (model in c('none','poisson')) for (d in 1:4) for (tie in c('sequence','support')) {
    a <- super_cluster2(tab, distance = d, tie_break = tie, indel_model = model, verbose = FALSE)
    b <- super_cluster2(tab, distance = d, tie_break = tie, indel_model = model, verbose = FALSE, use_kmer_filter = FALSE)
    expect_identical(canonical_members(a), canonical_members(b), info = paste(model,d,tie))
    expect_equal(sum(a$sum_counts), sum(tab$counts))
  }
})

test_that('bounded posting lists remain exact for unsorted native inputs', {
  set.seed(917)
  roots <- replicate(30, paste(sample(c('A','C','G','T'), 12, TRUE), collapse = ''))
  bc <- unique(c(roots, paste0(substr(roots,1,7),substr(roots,9,12)),
                 paste0(substr(roots,1,6),'A',substr(roots,7,12)),
                 '', 'N', 'NN', 'AAAN', 'AAA'))
  counts <- sample(c(1L,2L,5L,20L,100L,10000L),length(bc),TRUE)
  canon <- function(x) sort(paste(x$central_barcode, x$sum_counts,
    vapply(x$all_barcodes, function(s) paste(sort(s),collapse='|'), character(1))))
  for (model in c(FALSE,TRUE)) for (d in 2:4) {
    native <- function(filter) barbac:::barbac_cpp_centroid_cluster_optimized(
      bc, counts, d, 'lv', use_kmer_filter=filter,
      use_indel_model=model, verbose=FALSE)
    expect_identical(canon(native(TRUE)),canon(native(FALSE)))
  }
})

test_that('wide LV guard pruning retains near parents and floor-boundary assignments', {
  set.seed(20200910)
  roots <- replicate(120,paste(sample(c('A','C','G','T'),15,TRUE),collapse=''))
  variants <- unique(c(roots,paste0('A',substr(roots,2,15)),
    paste0('AA',substr(roots,3,15)),substr(roots,1,14),
    'GTAAAAAAAAAACTG','GTAAAAAAAAACTG'))
  tab <- data.frame(barcode=variants,counts=sample(c(1L,2L,4L,5L,9L,10L,99L,100L,299L,300L,10000L),length(variants),TRUE))
  for(model in c('none','poisson')) for(d in 2:4) for(ratio in c(1,20,20.01)) {
    indexed <- super_cluster2(tab,distance=d,method='lv',merge_ratio=ratio,
      indel_model=model,tie_break='support',verbose=FALSE)
    exhaustive <- super_cluster2(tab,distance=d,method='lv',merge_ratio=ratio,
      indel_model=model,tie_break='support',verbose=FALSE,use_kmer_filter=FALSE)
    expect_identical(canonical_members(indexed),canonical_members(exhaustive),info=paste(model,d,ratio))
  }
})

test_that('Poisson exception excludes nonrepeated gaps and permits plausible insertions', {
  tab <- data.frame(barcode = c('ACGTACGTACGT', 'ACTACGTACGT'), counts=c(1000L,65L))
  r <- super_cluster2(tab, distance=1, indel_model='poisson', verbose=FALSE)
  expect_equal(nrow(r),2)
  tab$barcode <- c('GTAAAAAAAAAACTG','GTAAAAAAAAAAACTG')
  r <- super_cluster2(tab, distance=1, indel_model='poisson',error_rate=0.02,verbose=FALSE)
  expect_equal(r$central_barcode,tab$barcode[1])
  expect_equal(r$sum_counts,1065L)
})
