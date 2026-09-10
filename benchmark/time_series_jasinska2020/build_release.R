#!/usr/bin/env Rscript
kit <- normalizePath('benchmark/time_series_jasinska2020')
dest <- file.path(kit,'generated/release');lib <- file.path(dest,'lib')
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
args <- commandArgs(TRUE)
if(!identical(args,'--record-existing')) {
  tar <- pkgbuild::build('.',dest_path=dest,vignettes=FALSE,manual=FALSE)
  log <- file.path(dest,'install.log')
  status <- system2(file.path(R.home('bin'),'R'),c('CMD','INSTALL','--preclean','--no-multiarch',
    '-l',shQuote(lib),shQuote(tar)),stdout=log,stderr=log)
  if(status!=0)stop('Release installation failed; see ',log)
}
log <- file.path(dest,'install.log')
stopifnot(file.exists(log))
compile <- grep(' -c clustering.cpp ',readLines(log,warn=FALSE),value=TRUE)
stopifnot(length(compile)==1,grepl(' -O[123s] ',compile),!grepl(' -O0 ',compile))
binary <- list.files(file.path(lib,'barbac/libs'),pattern='\\.(so|dll)$',full.names=TRUE,recursive=TRUE)
stopifnot(length(binary)==1)
sources <- c(list.files('R',pattern='\\.R$',full.names=TRUE),
  list.files('src',pattern='\\.(cpp|h)$',full.names=TRUE),'DESCRIPTION','NAMESPACE')
sha <- function(f)digest::digest(file=f,algo='sha256')
jsonlite::write_json(list(status='release_build',built_at_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
  compile_command=compile,binary_sha256=sha(binary),install_log_sha256=sha(log),
  source_hashes=as.list(setNames(vapply(sources,sha,character(1)),sources)),
  R=R.version.string,platform=R.version$platform),file.path(dest,'build.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Verified optimized release build in ',lib,'\n',sep='')
