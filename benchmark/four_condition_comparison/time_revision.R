# Repeat clustering in one loaded R session. Run revisions sequentially on an
# otherwise idle machine; this measures warm algorithm time, not process startup.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5L) stop('Expected library, data_root, output.csv, tie_break, repetitions')
.libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages(library(barbac))
rows <- list()
conditions <- c('random_substitutions','random_low_indels','anchored_substitutions','anchored_low_indels')
for (condition in conditions) {
  input <- file.path(args[[2]], condition, 'input.csv')
  for (method in c('hamming','lv')) {
    reference <- NULL
    for (iteration in seq_len(as.integer(args[[5]]))) {
      gc()
      elapsed <- system.time(result <- super_cluster2(input, method=method,
        tie_break=args[[4]], verbose=FALSE))[["elapsed"]]
      compact <- result[c('central_barcode','sum_counts')]
      if (is.null(reference)) reference <- compact else stopifnot(identical(reference,compact))
      row <- data.frame(condition=condition,method=method,tie_break=args[[4]],
        repeat_id=iteration,algorithm_s=elapsed,n_centroids=nrow(result),
        reads=sum(result$sum_counts),build_id=barbac:::barbac_build_id())
      rows[[length(rows)+1L]] <- row
      write.csv(do.call(rbind,rows),args[[3]],row.names=FALSE)
      cat(condition,method,iteration,elapsed,'seconds\n')
    }
  }
}
