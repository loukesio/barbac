#!/usr/bin/env Rscript
# Run from the repository root. CLI tools, app dependencies and Quarto are needed.
# New receipts are separate from the archived scientific validation records.
repo <- normalizePath('.')
out <- file.path(repo, 'benchmark', 'validation', 'release_2026-09-11')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
verify <- function() {
  pkgload::load_all(repo, quiet = TRUE)
  available <- barbac::check_barbac_tools(tools = c('fastqc', 'pear', 'minimap2', 'samtools'))
  stopifnot(all(available$available))
  results <- as.data.frame(testthat::test_dir('tests/testthat', reporter = 'summary', stop_on_failure = TRUE))
  stopifnot(sum(results$failed) == 0, !any(results$error), !any(results$skipped))
  rscript <- file.path(R.home('bin'), 'Rscript')
  stopifnot(system2(rscript, c('app/run.R', '--prepare-only')) == 0L)
  stopifnot(system2(rscript, 'app/tests/run_tests.R') == 0L)
  app <- jsonlite::read_json('app/.qa/engine-validation.json', simplifyVector = TRUE)
  stopifnot(app$status == 'passed', app$failures == 0, app$errors == 0, app$skipped == 0)
  sources <- c('DESCRIPTION', 'NAMESPACE', list.files('R', '[.]R$', full.names = TRUE),
    list.files('src', '[.](cpp|h)$', full.names = TRUE),
    list.files('tests/testthat', '[.]R$', full.names = TRUE),
    'app/app.R', 'app/run.R', list.files('app/R', '[.]R$', full.names = TRUE),
    'app/tests/test-engine.R', 'app/report.qmd', 'tools/verify_release.R')
  sha <- function(path) digest::digest(file = path, algo = 'sha256')
  list(status = 'passed', checked_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
    package_version = as.character(utils::packageVersion('barbac')),
    native_build = barbac:::barbac_build_id(),
    R = R.version.string, platform = R.version$platform,
    package_tests = list(cases = nrow(results), assertions = sum(results$passed),
      failures = sum(results$failed), errors = sum(results$error),
      skipped = sum(results$skipped), warnings = sum(results$warning)),
    app_tests = app, tools = available[c('tool', 'version')],
    integration = list(R1_only_without_PEAR = TRUE, mixed_single_paired = TRUE,
      expected_reads_per_fixture_sample = 144, expected_clusters_per_fixture_sample = 3,
      known_indel_sequences_and_counts = TRUE, fail_fast_commands = TRUE,
      native_membership_and_plot_parity = TRUE, Quarto_export = TRUE),
    source_hashes = as.list(setNames(vapply(sources, sha, character(1)), sources)),
    scope = 'Local correctness and integration verification using synthetic fixtures and real executables. Not a runtime benchmark or a validation of every library design.')
}
receipt <- tryCatch(verify(), error = function(e) {
  jsonlite::write_json(list(status = 'failed', error = conditionMessage(e)),
    file.path(out, 'validation.json'), pretty = TRUE, auto_unbox = TRUE)
  stop(e)
})
jsonlite::write_json(receipt, file.path(out, 'validation.json'), pretty = TRUE, auto_unbox = TRUE)
cat('Release verification passed:', file.path(out, 'validation.json'), '\n')
