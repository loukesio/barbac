#!/usr/bin/env Rscript
# Known barcodes at the read boundary, including indels and reverse orientation.
devtools::load_all(quiet=TRUE);barbac::use_barbac_env()
root <- normalizePath('benchmark/time_series_jasinska2020')
work <- tempfile('barbac-ecoli-extraction-');dir.create(work)
flank <- 'TATCTCGGTAGTGGGATACGACGATACCGAAGACA'
expected <- c('ACGACGACGACGACG','ACGACGACGACGAC','ACGACGACGACGACGA','ACGACGACGA')
seqs <- paste0(expected,flank)
seqs <- c(seqs,as.character(Biostrings::reverseComplement(Biostrings::DNAString(seqs[1]))),
          paste0(expected[1],'CATCTCGGTAGTGGGATACGACGATACCGAAGACA'))
names(seqs) <- paste0('read',seq_along(seqs))
fq <- file.path(work,'test.fastq')
writeLines(unlist(lapply(seq_along(seqs),function(i)c(paste0('@',names(seqs)[i]),seqs[i],'+',
  paste(rep('I',nchar(seqs[i])),collapse='')))),fq)
sam <- file.path(work,'test.sam');bam <- file.path(work,'test.bam')
stopifnot(system2('minimap2',c('-a','-x','sr','-k','9','-w','5','-m','10','-s','10','-n','1','--secondary=no',
  file.path(root,'reference/cassette.fasta'),fq),stdout=sam,stderr=file.path(work,'mapping.log'))==0L)
stopifnot(system2('samtools',c('sort','-o',bam,sam))==0L,system2('samtools',c('index',bam))==0L)
result <- barbac::barbac_xtr(bam,ref_name='Jasinska2020_barcode_cassette',start_pos=32,end_pos=40,
  flank_pattern='^([ACGT]{10,20})TATCTCGGTAG',include_read_ids=TRUE,
  output_file=file.path(work,'barcodes.csv'),verbose=FALSE)
calls <- readr::read_csv(result,show_col_types=FALSE)
calls <- calls[order(calls$read_id), ]
testthat::expect_identical(calls$read_id,paste0('read',1:5))
testthat::expect_identical(calls$barcode,c(expected,expected[1]))
testthat::expect_equal(calls$barcode_length,c(15,14,16,10,15))
jsonlite::write_json(list(status='passed',known_boundary_barcodes=5,
  insertion_deletion_and_reverse_orientation=TRUE,altered_flank_rejected=TRUE,
  reference_sha256=digest::digest(file=file.path(root,'reference/cassette.fasta'),algo='sha256')),
  file.path(root,'sources/extraction_validation.json'),pretty=TRUE,auto_unbox=TRUE)
unlink(work,recursive=TRUE)
cat('Known read-boundary, indel and orientation extraction checks passed.\n')
