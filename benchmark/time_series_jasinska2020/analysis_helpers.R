# Shared, directly testable calculations. No fitting to publication agreement.
diversity <- function(counts,total_reads) {
  stopifnot(all(is.finite(counts)),all(counts>=0),total_reads>0,sum(counts)<=total_reads)
  cts <- counts[counts>0]
  if(!length(cts)) return(c(richness=0,shannon_effective=NA,inverse_dominance=NA,
                          assigned_shannon_effective=NA,assigned_inverse_dominance=NA))
  p <- cts/total_reads
  q <- cts/sum(cts)
  c(richness=length(cts),shannon_effective=exp(-sum(p*log(p))),inverse_dominance=1/max(p),
    assigned_shannon_effective=exp(-sum(q*log(q))),assigned_inverse_dominance=1/max(q))
}
read_author_tables <- function(root) {
  f <- file.path(root,'sources/supplementary_tables.xlsx')
  s <- as.data.frame(suppressMessages(readxl::read_excel(f,sheet='Supplementary Table 1c',skip=1)))
  names(s) <- c('population','passage','generation','input_reads','quality_fraction',
                'extracted_fraction','raw_sequences','clustered_sequences')
  s <- s[!is.na(s$population), ]
  stopifnot(!anyDuplicated(s[c('population','passage')]))
  t <- as.data.frame(suppressMessages(readxl::read_excel(f,sheet='Supplementary Table 4b',skip=2,col_names=FALSE)))
  pop <- NA_character_; result <- list()
  for(i in seq_len(nrow(t))) {
    id <- t[i,1]
    if(is.na(id)) next
    if(grepl(' r[123]$',id)) {pop <- id; next}
    if(grepl('^[ACGT]{10,20}$',id)) result[[length(result)+1L]] <- data.frame(
      population=pop,barcode=id,published_mean=as.numeric(t[i,2]),
      published_final=as.numeric(t[i,3]),published_color=t[i,4])
  }
  top <- do.call(rbind,result)
  stopifnot(!anyDuplicated(top[c('population','barcode')]),all(top$published_final>=0),
            all(top$published_final<=1))
  list(samples=s,top=top)
}
collapse_counts <- function(frames) {
  z <- do.call(rbind,frames)
  out <- aggregate(z$counts,list(barcode=z$barcode),sum)
  names(out)[2] <- 'counts'
  out$barcode_length <- nchar(out$barcode)
  out[order(-out$counts,out$barcode), ]
}
