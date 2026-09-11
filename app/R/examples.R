studio_demo <- function() {
  set.seed(117)
  dna <- function(n) paste(sample(c('A','C','G','T'), n, TRUE), collapse='')
  barcodes <- replicate(18, dna(26))
  rows <- list()
  for (t in 0:7) {
    weights <- exp(seq(1.8, -.4, length.out=18) +
      sin(seq_len(18)*.8+t*.47)*.55 + c(rep(-.28,5),rep(.30,5),rep(.02,8))*t)
    counts <- round(weights/sum(weights)*24000)
    parent <- data.frame(sample=paste0('T',t),population='Example population',time=t*8,
                         barcode=barcodes,counts=counts)
    variants <- parent[1:8,]
    variants$barcode <- vapply(variants$barcode,function(z)
      paste0(substr(z,1,12),if(substr(z,13,13)=='A')'C' else 'A',substr(z,14,26)), character(1))
    variants$counts <- pmax(1, floor(variants$counts*.004))
    parent$counts[1:8] <- parent$counts[1:8]-variants$counts
    rows[[t+1]] <- rbind(parent,variants)
  }
  do.call(rbind,rows)
}

studio_fastq_example <- function(directory) {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  set.seed(20260910)
  dna <- function(n) paste(sample(c('A','C','G','T'),n,TRUE),collapse='')
  prefix <- dna(170); suffix <- dna(170); truth <- replicate(3,dna(26))
  bc <- rep(truth,c(80,40,20))
  molecules <- paste0(prefix,bc,suffix)
  writeLines(c('>example_cassette',paste0(prefix,strrep('N',26),suffix)),file.path(directory,'reference.fasta'))
  for (mate in 1:2) {
    seqs <- if(mate==1)substr(molecules,1,250) else
      substr(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(molecules))),1,250)
    con <- gzfile(file.path(directory,paste0('example_R',mate,'.fastq.gz')),'wt')
    writeLines(unlist(lapply(seq_along(seqs),function(i)c(paste0('@example_',i,'/',mate),seqs[i],'+',strrep('I',250)))),con)
    close(con)
  }
  # Single-end example spans the entire cassette.
  con <- gzfile(file.path(directory,'example_single.fastq.gz'),'wt')
  writeLines(unlist(lapply(seq_along(molecules),function(i)c(paste0('@single_',i),molecules[i],'+',strrep('I',nchar(molecules[i]))))),con);close(con)
  readr::write_csv(data.frame(barcode=truth,counts=c(80,40,20)),file.path(directory,'expected_barcodes.csv'))
  writeLines(c('Synthetic barbac Studio extraction example; 140 molecules, 3 known barcodes.',
    'Paired route: upload example_R1.fastq.gz and example_R2.fastq.gz.',
    'Single-end route: upload example_single.fastq.gz as R1; leave R2 empty.',
    'Reference: reference.fasta. Barcode start: 171. Barcode end: 196.',
    paste0('Left flank: ',substr(prefix,159,170)),paste0('Right flank: ',substr(suffix,1,12)),
    'Minimum length: 24. Maximum length: 28. Choose exact-flank extraction.',
    'These are synthetic sequencing reads, not a biological measurement.'),file.path(directory,'README.txt'))
  list(start=171,end=196,left=substr(prefix,159,170),right=substr(suffix,1,12),min_length=24,max_length=28,mode='flanks')
}
