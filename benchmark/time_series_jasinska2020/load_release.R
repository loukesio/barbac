# Use an isolated installed package: devtools::load_all() defaults to -O0.
load_release <- function(kit) {
  source(file.path(kit,'provenance_helpers.R'),local=TRUE)
  lib <- file.path(kit,'generated/release/lib')
  metadata <- file.path(kit,'generated/release/build.json')
  if(!file.exists(metadata))stop('Build the release package first: Rscript ',file.path(kit,'build_release.R'))
  build <- jsonlite::read_json(metadata,simplifyVector=TRUE)
  sha <- function(f)digest::digest(file=f,algo='sha256')
  assert_hash_map(build$source_hashes,file.path(kit,'../..'))
  .libPaths(c(lib,.libPaths()))
  library(barbac,lib.loc=lib,character.only=FALSE)
  binary <- getLoadedDLLs()[['barbac']][['path']]
  stopifnot(startsWith(normalizePath(binary),normalizePath(lib)),sha(binary)==build$binary_sha256)
  build
}
