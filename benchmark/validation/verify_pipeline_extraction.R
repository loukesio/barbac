#!/usr/bin/env Rscript
# Run from the repository root. Uses real CLI tools and independently known
# synthetic reads; no sequencing downloads or historical results are changed.
repo <- normalizePath('.')
out <- file.path(repo, 'benchmark/validation/pipeline_extraction')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(list(status = 'running'), file.path(out, 'validation.json'), auto_unbox = TRUE)
work <- tempfile('barbac-pipeline-verification-')
dir.create(work)

verify <- function() {
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  pkgload::load_all(repo, quiet = TRUE, export_all = FALSE)
  barbac::use_barbac_env()
  tools <- barbac::check_barbac_tools()
  stopifnot(all(tools$available))
  set.seed(20260910)
  dna <- function(n) paste(sample(c('A', 'C', 'G', 'T'), n, replace = TRUE), collapse = '')
  prefix <- dna(170)
  suffix <- dna(170)
  truth <- setNames(replicate(3, dna(26)), c('A', 'B', 'C'))
  a <- truth[['A']]
  variants <- c(insertion = paste0(substr(a, 1, 12), 'A', substr(a, 13, 26)),
    deletion = paste0(substr(a, 1, 12), substr(a, 14, 26)),
    substitution = paste0(substr(a, 1, 12), if (substr(a, 13, 13) == 'A') 'C' else 'A', substr(a, 14, 26)))
  truth_distance <- stringdist::stringdistmatrix(truth, truth, method = 'lv')
  stopifnot(all(truth_distance[upper.tri(truth_distance)] > 6))
  reference <- file.path(work, 'cassette.fasta')
  writeLines(c('>verification_cassette', paste0(prefix, strrep('N', 26), suffix)), reference)
  expected <- list(t0 = setNames(c(80, 40, 20), truth),
    t1 = setNames(c(40, 80, 20, 2, 2, 1), c(truth, variants)))
  expected_unmapped <- c(t0 = 0, t1 = 5)
  input_rows <- list()
  expected_reads <- list()
  for (sample in names(expected)) {
    bc <- rep(names(expected[[sample]]), expected[[sample]])
    sequences <- c(paste0(prefix, bc, suffix), rep(dna(366), expected_unmapped[[sample]]))
    read_ids <- paste0(sample, '_', seq_along(sequences))
    expected_reads[[sample]] <- setNames(bc, read_ids[seq_along(bc)])
    write_fastq <- function(path, seqs, mate) {
      records <- unlist(lapply(seq_along(seqs), function(i) c(
        paste0('@', read_ids[i], '/', mate), seqs[i], '+', strrep('I', nchar(seqs[i])))))
      writeLines(records, path)
    }
    r1 <- file.path(work, paste0(sample, '_R1.fastq'))
    r2 <- file.path(work, paste0(sample, '_R2.fastq'))
    write_fastq(r1, substr(sequences, 1, 250), 1)
    reverse <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequences)))
    write_fastq(r2, substr(reverse, 1, 250), 2)
    input_rows[[sample]] <- data.frame(sample = sample, R1 = r1, R2 = r2)
  }
  samples <- do.call(rbind, input_rows)
  pipe <- barbac::run_cli_pipeline(samples, reference, file.path(work, 'results'), verbose = FALSE)
  stopifnot(nrow(pipe$stats) == 2)
  count_rows <- list()
  extracted <- list()
  for (sample in names(expected)) {
    bam <- file.path(pipe$bam_dir, paste0(sample, '_ANC.assembled_sorted.bam'))
    stopifnot(file.exists(bam), file.exists(paste0(bam, '.bai')))
    mapping <- pipe$stats[pipe$stats$sample == paste0(sample, '_ANC.assembled'), ]
    stopifnot(nrow(mapping) == 1, mapping$mapped == sum(expected[[sample]]),
      mapping$unmapped == expected_unmapped[[sample]])
    assembled <- file.path(pipe$merged_dir, paste0(sample, '_ANC.assembled.fastq'))
    stopifnot(length(readLines(assembled)) / 4 == sum(expected[[sample]]) + expected_unmapped[[sample]])
    for (mate in 1:2) {
      qc <- file.path(pipe$fastqc_dir, paste0(sample, '_R', mate, '_fastqc.zip'))
      stopifnot(file.exists(qc))
      con <- unz(qc, paste0(sample, '_R', mate, '_fastqc/fastqc_data.txt'))
      qc_text <- readLines(con); close(con)
      total <- as.integer(sub('Total Sequences\t', '', grep('^Total Sequences\t', qc_text, value = TRUE)))
      stopifnot(total == sum(expected[[sample]]) + expected_unmapped[[sample]])
    }
    pattern <- paste0(substr(prefix, 159, 170), '([ACGT]{24,28})', substr(suffix, 1, 12))
    csv <- barbac::barbac_xtr(bam, 'verification_cassette', 171, 196,
      output_file = file.path(work, paste0(sample, '_barcodes.csv')),
      flank_pattern = pattern, yield_size = 17, verbose = FALSE)
    x <- read.csv(csv)
    stopifnot(setequal(x$barcode, names(expected[[sample]])),
      identical(as.numeric(x$counts[match(names(expected[[sample]]), x$barcode)]),
                as.numeric(expected[[sample]])), all(x$barcode_length == nchar(x$barcode)))
    read_csv <- barbac::barbac_xtr(bam, 'verification_cassette', 171, 196,
      output_file = file.path(work, paste0(sample, '_read_ids.csv')),
      flank_pattern = pattern, include_read_ids = TRUE, yield_size = 13, verbose = FALSE)
    reads <- read.csv(read_csv)
    # PEAR retains the /1 suffix from each read-pair identifier.
    ids <- sub('/1$', '', reads$read_id)
    stopifnot(!anyDuplicated(ids), setequal(ids, names(expected_reads[[sample]])),
      identical(reads$barcode, unname(expected_reads[[sample]][ids])))
    extracted[[sample]] <- x
    count_rows[[sample]] <- data.frame(sample = sample, barcode = names(expected[[sample]]),
      length = nchar(names(expected[[sample]])), expected = as.numeric(expected[[sample]]),
      observed = x$counts[match(names(expected[[sample]]), x$barcode)])
    if (sample == 't0') {
      fixed <- barbac::barbac_xtr(bam, 'verification_cassette', 171, 196,
        output_file = file.path(work, 'fixed.csv'), verbose = FALSE)
      fixed <- read.csv(fixed)
      stopifnot(setequal(fixed$barcode, names(expected[[sample]])),
        identical(as.numeric(fixed$counts[match(names(expected[[sample]]), fixed$barcode)]),
                  as.numeric(expected[[sample]])))
    }
  }
  stopifnot(file.exists(file.path(pipe$output_dir, 'multiqc/multiqc_report.html')),
    !any(grepl('failed with exit code|Command failed|FASTQ file not found', readLines(pipe$log_file))))
  pooled <- do.call(rbind, extracted)
  pooled <- stats::aggregate(counts ~ barcode, pooled, sum)
  clusters <- barbac::super_cluster2(pooled, method = 'lv', distance = 3,
    tie_break = 'support', merge_ratio = 20, error_rate = 0.005, indel_model = 'poisson', verbose = FALSE)
  stopifnot(nrow(clusters) == 3, setequal(clusters$central_barcode, truth),
    identical(as.numeric(clusters$sum_counts[match(truth, clusters$central_barcode)]), c(125, 120, 40)))
  membership <- do.call(rbind, lapply(seq_len(nrow(clusters)), function(i)
    data.frame(observed_barcode = clusters$all_barcodes[[i]], barcode = clusters$central_barcode[i])))
  stopifnot(!anyDuplicated(membership$observed_barcode))
  counts <- do.call(rbind, lapply(names(extracted), function(sample) {
    x <- extracted[[sample]]
    data.frame(barcode = membership$barcode[match(x$barcode, membership$observed_barcode)],
      time = if (sample == 't0') 0 else 1, counts = x$counts)
  }))
  counts <- stats::aggregate(counts ~ barcode + time, counts, sum)
  target <- data.frame(barcode = rep(unname(truth), 2), time = rep(0:1, each = 3),
    expected = c(80, 40, 20, 45, 80, 20))
  target <- merge(target, counts, by = c('barcode', 'time'))
  stopifnot(nrow(target) == 6, all(target$expected == target$counts))
  qc <- barbac::barbac_xtr.stats(extracted$t1, c(25, 27), panel_labels = TRUE,
    return_details = TRUE, verbose = FALSE)
  stopifnot(sum(qc$length_summary$reads) == 145, sum(qc$length_summary$barcodes) == 6)
  cluster_qc <- barbac::cluster_stats(clusters, verbose = FALSE)
  stopifnot(cluster_qc$total_reads == 285, cluster_qc$n_clusters == 3)
  plot <- barbac::barbac_ts_area(counts, min_total_count = 0, fill_missing = 'zero',
    palette = 'alger', theme = barbac::theme_barbac(family = 'sans'))
  stopifnot(all(abs(tapply(plot$data$.freq, plot$data$time, sum) - 1) < 1e-12),
    length(unique(ggplot2::ggplot_build(plot)$data[[1]]$group)) == 3)
  ggplot2::ggsave(file.path(out, 'extraction_diagnostics.pdf'), qc$plot, width = 10, height = 7)
  ggplot2::ggsave(file.path(out, 'barcode_time_series.pdf'), plot, width = 7, height = 4)
  ggplot2::ggsave(file.path(out, 'mapping_qc.pdf'), barbac::plot_bam_stats(pipe$stats), width = 7, height = 4)
  readr::write_csv(do.call(rbind, count_rows), file.path(out, 'barcode_counts.csv'))
  readr::write_csv(target, file.path(out, 'lineage_counts.csv'))
  readr::write_csv(pipe$stats, file.path(out, 'bam_summary.csv'))
  cat('Real FASTQ -> PEAR -> minimap2/samtools -> extraction -> clustering -> plots passed.\n')
  # The unit suite also exercises internal helpers through the development loader.
  pkgload::load_all(repo, quiet = TRUE)
  results <- as.data.frame(testthat::test_dir(file.path(repo, 'tests/testthat'), reporter = 'summary', stop_on_failure = TRUE))
  stopifnot(sum(results$failed) == 0, !any(results$error), !any(results$skipped))
  sources <- c('DESCRIPTION', 'NAMESPACE', list.files('R', full.names = TRUE),
    list.files('src', pattern = '[.](cpp|h)$', full.names = TRUE),
    list.files('tests/testthat', pattern = '[.]R$', full.names = TRUE),
    'benchmark/validation/verify_pipeline_extraction.R')
  sha <- function(f) digest::digest(file = f, algo = 'sha256')
  artifacts <- c('barcode_counts.csv', 'lineage_counts.csv', 'bam_summary.csv',
    'extraction_diagnostics.pdf', 'barcode_time_series.pdf', 'mapping_qc.pdf')
  list(status = 'passed', checked_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
    base_revision = system2('git', c('rev-parse', 'HEAD'), stdout = TRUE),
    R = R.version.string, platform = R.version$platform, tools = tools[c('tool', 'version')],
    fixture = list(seed = 20260910, samples = 2, FASTQs = 4, input_pairs = 290,
      merged_reads = 290, mapped_reads = 285, unmapped_reads = 5, extracted_reads = 285,
      fixed_coordinate_exact_counts = TRUE, flank_exact_counts_and_read_ids = TRUE,
      observed_lengths = c(25, 26, 27), expected_centroids = 3, all_time_series_cells_match = TRUE),
    checks = list(FastQC_counts = TRUE, PEAR_counts = TRUE, BAM_counts = TRUE,
      MultiQC_report = TRUE, extraction_diagnostics = TRUE, clustering_conserves_reads = TRUE,
      native_palette_plot = TRUE),
    unit_tests = list(cases = nrow(results), assertions_passed = sum(results$passed),
      failures = sum(results$failed), errors = sum(results$error), skipped = sum(results$skipped),
      warnings = sum(results$warning)),
    source_hashes = as.list(setNames(vapply(sources, sha, character(1)), sources)),
    artifact_hashes = as.list(setNames(vapply(file.path(out, artifacts), sha, character(1)), artifacts)),
    scope = 'Small paired-end synthetic integration check using real executables; not an accuracy benchmark or verification of every input design.')
}
receipt <- tryCatch(verify(), error = function(e) {
  jsonlite::write_json(list(status = 'failed', error = conditionMessage(e)),
    file.path(out, 'validation.json'), pretty = TRUE, auto_unbox = TRUE)
  stop(e)
})
jsonlite::write_json(receipt, file.path(out, 'validation.json'), pretty = TRUE, auto_unbox = TRUE)
cat('Verification receipt: ', file.path(out, 'validation.json'), '\n', sep = '')
