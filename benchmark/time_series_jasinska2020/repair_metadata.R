#!/usr/bin/env Rscript
# One-time migration of the September 2026 local receipts. No numerical rerun:
# verify the original signed cache inputs, preserve timings, then retain JSON keys.
kit <- normalizePath('benchmark/time_series_jasinska2020')
source(file.path(kit,'provenance_helpers.R'))
backup <- file.path(kit,'generated/provenance_schema1')
state <- readRDS(file.path(backup,'state.rds'))
sha <- function(f)digest::digest(file=f,algo='sha256')
build <- jsonlite::read_json(file.path(backup,'build.json'),simplifyVector=TRUE)
sources <- c(list.files('R',pattern='\\.R$',full.names=TRUE),
  list.files('src',pattern='\\.(cpp|h)$',full.names=TRUE),'DESCRIPTION','NAMESPACE')
current <- setNames(vapply(sources,sha,character(1)),sources)
stopifnot(identical(unname(current),unname(build$source_hashes)))
build$source_hashes <- as.list(current)
assert_hash_map(build$source_hashes,'.')
binary <- list.files(file.path(kit,'generated/release/lib/barbac/libs'),pattern='\\.(so|dll)$',recursive=TRUE,full.names=TRUE)
stopifnot(length(binary)==1,sha(binary)==build$binary_sha256)
jsonlite::write_json(build,file.path(kit,'generated/release/build.json'),pretty=TRUE,auto_unbox=TRUE)
wells <- c('A3','B3','C3','A1','B1','C1');checks <- list()
for(w in wells) {
  original <- readRDS(file.path(backup,paste0(w,'.rds')))
  assert_hash_map(as.list(original$input_hashes),file.path(kit,'generated/clustering',w))
  expected <- digest::digest(list(original$input_hashes,build$binary_sha256,
    state$hashes[['benchmark_release.R']]),algo='sha256')
  stopifnot(identical(original$signature,expected),identical(original$status,'passed'))
  migrated <- original
  migrated$signature <- release_signature(original$input_hashes,build)
  test_old <- original;test_new <- migrated
  test_old$signature <- test_new$signature <- NULL
  stopifnot(identical(test_old,test_new))
  saveRDS(migrated,file.path(kit,'generated/release/checks',paste0(w,'.rds')))
  checks[[w]] <- list(original_signature=original$signature,new_signature=migrated$signature,
    all_non_signature_fields_identical=TRUE)
}
stopifnot(sha(file.path(backup,'release_clustering_times.csv'))==
  sha(file.path(kit,'results/release_clustering_times.csv')))
jsonlite::write_json(list(status='passed',reason='Retain filename keys in JSON checksum maps; strengthen numerical cache identity',
  numerical_package_sources_unchanged=TRUE,all_six_original_cache_signatures_verified=TRUE,
  timings_unchanged=TRUE,timing_csv_sha256=sha(file.path(kit,'results/release_clustering_times.csv')),
  before_script_hashes=as.list(state$hashes),
  after_script_hashes=as.list(setNames(vapply(names(state$hashes),function(f)sha(file.path(kit,f)),character(1)),names(state$hashes))),
  settings=release_settings(),populations=checks),
  file.path(kit,'results/metadata_schema_migration.json'),pretty=TRUE,auto_unbox=TRUE)
cat('All original cache signatures verified; every numerical value and package source preserved.\n')
