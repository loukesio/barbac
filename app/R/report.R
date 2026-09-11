studio_report <- function(result, destination, template, palette = 'alger') {
  quarto <- Sys.which('quarto')
  if (!nzchar(quarto)) stop('Install the Quarto CLI to export an HTML report.',call.=FALSE)
  folder <- tempfile('studio-report-');dir.create(folder)
  on.exit(unlink(folder,recursive=TRUE),add=TRUE)
  saveRDS(result,file.path(folder,'analysis.rds'))
  file.copy(template,file.path(folder,'report.qmd'))
  populations <- unique(result$time_series$population)
  plots <- list()
  for(i in seq_along(populations)) {
    p <- tryCatch(studio_area(result,populations[i],palette),error=function(e)NULL)
    if(!is.null(p)) {
      name <- paste0('lineages-',i,'.png')
      ggplot2::ggsave(file.path(folder,name),p,width=10,height=4.6,dpi=150,bg='white')
      plots[[length(plots)+1]] <- list(population=populations[i],file=name)
    }
  }
  saveRDS(plots,file.path(folder,'plots.rds'))
  log <- file.path(folder,'render.log')
  code <- system2(quarto,c('render',shQuote(file.path(folder,'report.qmd')),'--to','html','--quiet'),stdout=log,stderr=log)
  if(code!=0 || !file.exists(file.path(folder,'report.html')))
    stop('Report rendering failed: ',paste(tail(readLines(log,warn=FALSE),6),collapse=' '),call.=FALSE)
  file.copy(file.path(folder,'report.html'),destination,overwrite=TRUE)
  invisible(destination)
}
