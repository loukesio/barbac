# this is useful for counting reads in bin of the bam files

#https://combine-australia.github.io/2017-05-19-bioconductor-melbourne/AdvGRanges_Rtracklayer_Rsamtools.html

which <- GRanges(seqnames = c("Reference_barcodes"),
                 ranges = IRanges(c(1),c(174)))

what <- c("rname", "strand", "pos", "qwidth", "seq")

flag = scanBamFlag(isDuplicate = FALSE, isProperPair = TRUE, 
                   isPaired = TRUE)

param <- ScanBamParam(which = which, what = what, flag=flag)

reads = scanBam("test_4_seqs_sorted.bam", param = param)


names(reads)

param <- ScanBamParam(which = which, what = what)
param
test.bam

test.bam <- scanBam("test_4_seqs_sorted.bam", param = param)

a=test.bam$`Reference_barcodes:1-174`$seq

gr <- GRanges(seqnames = 'Reference_barcodes', ranges = IRanges(start = c(54, 78), width = 1))

library(BSgenome)
getSeq(test.bam, GRanges("Reference_barcodes:54-78"))

library("chromstaR")
reads <- readBamFileAsGRanges("test_4_seqs_sorted.bam", chromosomes='Reference_barcodes', pairedEndReads=FALSE,
                              min.mapq=10, remove.duplicate.reads=TRUE)

reads

getSeq(reads, GRanges("Reference_barcodes:54-78"))