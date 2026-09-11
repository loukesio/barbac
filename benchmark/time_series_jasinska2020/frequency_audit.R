# Keep both denominators visible; the author's final-frequency normalization
# cannot be reconciled completely from the available supplementary tables.
compare_frequency_denominators <- function(top,samples) {
  end <- samples[samples$passage==ave(samples$passage,samples$population,FUN=max),]
  stopifnot(!anyDuplicated(end$population),setequal(top$population,end$population))
  frac <- end$extracted_reads/end$input_reads_barbac
  stopifnot(all(is.finite(frac)),all(frac>0 & frac<=1))
  top$barbac_final_extracted <- top$barbac_final/frac[match(top$population,end$population)]
  stopifnot(all(top$barbac_final_extracted>=0 & top$barbac_final_extracted<=1))
  summary <- audit <- list()
  for(pop in unique(top$population)) {
    x <- top[top$population==pop,];e <- end[end$population==pop,]
    summary[[pop]] <- data.frame(Population=pop,Published_barcodes=nrow(x),
      Spearman=cor(x$published_final,x$barbac_final,method='spearman'),
      Mean_absolute_difference_input_pp=100*mean(abs(x$barbac_final-x$published_final)),
      Mean_absolute_difference_extracted_pp=100*mean(abs(x$barbac_final_extracted-x$published_final)),
      Maximum_absolute_difference_input_pp=100*max(abs(x$barbac_final-x$published_final)),
      Maximum_absolute_difference_extracted_pp=100*max(abs(x$barbac_final_extracted-x$published_final)))
    audit[[pop]] <- data.frame(Population=pop,Passage=e$passage,
      Published_top20_sum=sum(x$published_final),Published_extracted_fraction=e$extracted_fraction,
      Barbac_extracted_fraction=e$extracted_reads/e$input_reads_barbac,
      Incompatible_with_all_input=sum(x$published_final)>e$extracted_fraction+0.00005+1e-12)
  }
  list(top=top,summary=do.call(rbind,summary),audit=do.call(rbind,audit))
}
