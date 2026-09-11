#!/usr/bin/env Rscript
# Build reference documentation and copy only the synthetic Studio presentation.
# No runtime library, uploads, local reports or internal notes enter the site.
repo <- normalizePath('.')
source_dir <- tempfile('barbac-public-docs-')
dir.create(source_dir)
build_public_site <- function() {
  on.exit(unlink(source_dir, recursive = TRUE), add = TRUE)
  # Build from an allowlist of versioned package sources. pkgdown scans root
  # Markdown files, so .gitignore alone cannot exclude private working notes.
  tracked <- system2('git', c('ls-files'), stdout = TRUE)
  roots <- c('DESCRIPTION', 'NAMESPACE', 'README.md', 'NEWS.md', 'LICENSE.md',
             '.Rbuildignore', '_pkgdown.yml')
  sources <- tracked[tracked %in% roots | grepl('^(R|src|man|inst|vignettes)/', tracked)]
  for (path in sources) {
    destination <- file.path(source_dir, path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(repo, path), destination, overwrite = TRUE))
  }
  # pkgdown expects README illustrations under man/figures. Adapt only this
  # temporary copy; GitHub keeps the app's canonical presentation assets.
  preview <- file.path(source_dir, 'man', 'figures', 'studio-preview.gif')
  stopifnot(file.copy(file.path(repo, 'app/media/studio-preview.gif'), preview))
  readme <- file.path(source_dir, 'README.md')
  text <- readLines(readme, warn = FALSE)
  # The website template already displays the package logo beside its title.
  text <- text[!grepl('^<img align="right" src="man/figures/logo.png"', text)]
  writeLines(gsub('app/media/studio-preview.gif', 'man/figures/studio-preview.gif',
                  text, fixed = TRUE), readme)
  pkgdown::build_site_github_pages(pkg = source_dir, dest_dir = file.path(repo, 'docs'),
    new_process = FALSE, install = FALSE)
}
build_public_site()
media <- c('studio-preview.gif', 'studio-walkthrough.mp4', 'studio-fastq.mp4',
           'studio-walkthrough.jpg', 'studio-fastq.jpg', 'watch.html', 'recording.json')
target <- file.path('docs', 'app', 'media')
dir.create(target, recursive = TRUE, showWarnings = FALSE)
stopifnot(all(file.copy(file.path('app', 'media', media), target, overwrite = TRUE)))
# These guides and scientific scripts are viewed as source on GitHub. Preserve
# relative media URLs above so videos play directly from the documentation site.
index <- file.path('docs', 'index.html')
html <- readLines(index, warn = FALSE)
matches <- gregexpr('href="[^"]+"', html)
links <- regmatches(html, matches)
regmatches(html, matches) <- lapply(links, function(line) vapply(line, function(link) {
  href <- substring(link, 7, nchar(link) - 1)
  if (!grepl('^((documentation|benchmark|[.]github|inst)/|app/README[.])', href)) return(link)
  path <- sub('#.*$', '', href)
  fragment <- substring(href, nchar(path) + 1)
  # pkgdown rewrites Markdown links to .html even for source-only guides.
  markdown <- sub('[.]html$', '.md', path)
  source <- if (file.exists(file.path(repo, path))) path else markdown
  if (!file.exists(file.path(repo, source))) stop('Unresolved documentation source: ', href)
  paste0('href="https://github.com/loukesio/barbac/blob/main/', source, fragment, '"')
}, character(1)))
writeLines(html, index)
