# JSON objects must retain filename keys; unnamed arrays cannot validate files.
assert_hash_map <- function(hashes,root) {
  keys <- names(hashes)
  stopifnot(length(hashes)>0,length(keys)==length(hashes),
    !anyNA(keys),all(nzchar(keys)),!anyDuplicated(keys))
  for(f in keys) {
    expected <- hashes[[f]]
    stopifnot(is.character(expected),length(expected)==1,
      grepl('^[a-f0-9]{64}$',expected),file.exists(file.path(root,f)),
      identical(digest::digest(file=file.path(root,f),algo='sha256'),expected))
  }
  invisible(TRUE)
}
release_settings <- function()list(method='lv',distance=3,tie_break='support',
  merge_ratio=20,error_rate=.005,indel_model='poisson',use_design=FALSE)
release_signature <- function(hashes,build)digest::digest(list(input=hashes,
  binary=build$binary_sha256,package_sources=build$source_hashes,
  settings=release_settings()),algo='sha256')
