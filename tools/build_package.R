#!/usr/bin/env Rscript
# Build the package from versioned package files, without traversing local raw data.
repo <- normalizePath('.')
build_package <- function() {
  source <- tempfile('barbac-package-')
  dir.create(source)
  on.exit(unlink(source, recursive = TRUE), add = TRUE)
  paths <- system2('git', 'ls-files', stdout = TRUE)
  roots <- c('DESCRIPTION', 'NAMESPACE', 'README.md', 'NEWS.md', 'LICENSE.md', '.Rbuildignore')
  paths <- paths[paths %in% roots | grepl('^(R|src|man|inst|tests|vignettes)/', paths)]
  for (path in paths) {
    destination <- file.path(source, path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(repo, path), destination))
  }
  destination <- file.path(repo, 'app', '.qa', 'release')
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  pkgbuild::build(source, dest_path = destination, vignettes = TRUE, manual = FALSE)
}
build_package()
