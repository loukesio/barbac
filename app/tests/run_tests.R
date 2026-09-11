#!/usr/bin/env Rscript
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
app<-dirname(dirname(normalizePath(script)));repo<-dirname(app)
setwd(repo)
results<-as.data.frame(testthat::test_file('app/tests/test-engine.R',reporter='summary',stop_on_failure=TRUE))
dir.create('app/.qa',showWarnings=FALSE)
jsonlite::write_json(list(status='passed',cases=nrow(results),assertions=sum(results$passed),
  failures=sum(results$failed),errors=sum(results$error),skipped=sum(results$skipped)),
  'app/.qa/engine-validation.json',pretty=TRUE,auto_unbox=TRUE)
