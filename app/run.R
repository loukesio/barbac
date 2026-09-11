#!/usr/bin/env Rscript
# Launch from any working directory. Keep a release build separate from the
# user's installed barbac and from the archived scientific analyses.
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])
here <- dirname(normalizePath(script))
repo <- dirname(here)
required <- c('shiny', 'bslib', 'DT', 'future', 'promises', 'jsonlite',
              'digest', 'ggiraph', 'zip')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Install app dependencies: install.packages(c(',
                         paste(sprintf('"%s"', missing), collapse = ', '), '))')
if (packageVersion('shiny') < '1.8.1') stop('Shiny >= 1.8.1 is required.')
runtime <- file.path(here, '.runtime')
lib <- file.path(runtime, 'library')
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
sources <- c(file.path(repo, c('DESCRIPTION', 'NAMESPACE')),
  list.files(file.path(repo, 'R'), '[.]R$', full.names = TRUE),
  list.files(file.path(repo, 'src'), '[.](cpp|h)$', full.names = TRUE))
signature <- digest::digest(vapply(sources, digest::digest, character(1),
                                  algo = 'sha256', file = TRUE), algo = 'sha256')
stamp <- file.path(runtime, 'source.sha256')
if (!file.exists(stamp) || readLines(stamp, warn = FALSE)[1] != signature ||
    !file.exists(file.path(lib, 'barbac', 'DESCRIPTION'))) {
  message('Preparing the current barbac release for Studio…')
  log <- file.path(runtime, 'install.log')
  status <- system2(file.path(R.home('bin'), 'R'),
    c('CMD', 'INSTALL', '--preclean', '--no-multiarch', '--no-docs',
      paste0('--library=', shQuote(lib)), shQuote(repo)), stdout = log, stderr = log)
  if (status != 0L) stop(paste(readLines(log, warn = FALSE), collapse = '\n'))
  writeLines(signature, stamp)
}
.libPaths(c(lib, .libPaths()))
Sys.setenv(BARBAC_STUDIO_SOURCE_SHA256 = signature)
if ('--prepare-only' %in% commandArgs(TRUE)) {
  message('Studio release library is ready: ', lib)
  quit(status = 0)
}
port <- as.integer(Sys.getenv('BARBAC_STUDIO_PORT', '3838'))
if(Sys.getenv('BARBAC_STUDIO_TRACE')=='1')options(shiny.trace=TRUE)
shiny::runApp(here, host = Sys.getenv('BARBAC_STUDIO_HOST', '127.0.0.1'),
              port = port, launch.browser = interactive())
