################################
# libraries 
################################

library(here)         # set the working directory 
library(tidyverse)    # best package ever
library(Biostrings)   # read sequence data; be careful of the overlap with tidyverse

###############################
# set the working directory
###############################

setwd(here("data"))

###############################
# read data 
###############################

#$ the data contains 4 sequences that have been mapped to the reference

# extract the sequnce after using str_extract
(test <- readDNAStringSet("covert_bam_fasta.fasta") %>%   
  reverseComplement() %>%        # I do not understand why we need the reverse complement 
  as.data.frame() %>%
  as_tibble() %>% 
 dplyr::rename(seq=x) %>% 
  mutate(seq=str_extract(string = .$seq, pattern = "(?<=CTACG).*(?=CAGTC)"))
)


setwd(here("Volumes","home","test_bam"))
setwd("/Volumes/home/test_bam")

(test <- readDNAStringSet("20_3_shaking.fasta") %>%   
    reverseComplement() %>%        # I do not understand why we need the reverse complement 
    as.data.frame() %>%
    as_tibble() %>% 
    dplyr::rename(seq=x) %>% 
    mutate(seq=str_extract(string = .$seq, pattern = "(?<=CTACG).*(?=CAGTC)"))
)



test2 <- readDNAStringSet("covert_bam_fasta.fasta") %>%   
  reverseComplement()
test2

trimLRPatterns(Lpattern = "CTACG", Rpattern = "CAGTC", subject = test2[[1]],
               max.Lmismatch = 0, max.Rmismatch = 0)



trim_outside_patterns <- function(sequ, Lpattern, Rpattern){
  matched_view <- matchLRPatterns(Lpattern, Rpattern,
                                  subject = sequ,
                                  max.gaplength = 30,
                                  max.Lmismatch = 1,
                                  max.Rmismatch = 1,
                                  with.Lindels = TRUE,
                                  with.Rindels = TRUE)
  # select longest match
  
  matched_view[[which.max(width(matched_view))]]
}

trimmed <- endoapply(test2, trim_outside_patterns, "CTACG", "CAGTC")


trimLRPatterns("CTACG", "CAGTC",
               subject = test2[[1]],
               max.Lmismatch = 0,
               max.Rmismatch = 0,
               with.Lindels = TRUE,
               with.Rindels = TRUE) 


check <-c()
for(i in 1:5){
  a=matchLRPatterns("CTACG", "CAGTC", i,max.Rmismatch=0, test2)
  bind <- cbind(a,i)
  check<-rbind(check,bind)
}


matchLRPatterns("CTACG", "CAGTC",35, test2[[1]])

agrep("laysy", c("1 lazy", "1", "1 LAZY"), max.distance = 2)


#https://stat.ethz.ch/R-manual/R-devel/library/base/html/agrep.html



matchLRPatterns(Lpattern, Rpattern, 500, subject) # 1 match

Lpattern <- "ATTGCGCGA"
Rpattern <- "CGAAAATTTA"


dna_string <- DNAString("AAAAAAAANNNNNNNNNNNNNNNNNNNNNNNNNCCCCC")
dna_string

# I am writing this for loop in here to run over the differnet lengths and then 
# idea is that at output I ll erase the first 5 and last 5 characters 
# after I do this i need to keep the first non NA result in the for loop and I 
# need to find a way how to do it 
check <-c()
for(i in 6:30){
  a=matchLRPatterns("AAAAA", "CCCCC", i,max.Rmismatch=0, dna_string)
  bind <- cbind(a,i)
  check<-rbind(check,bind)
}

a=matchLRPatterns("AAAAA", "CCCCC", 25,max.Rmismatch=0, dna_string)
a
str(a)
a@subject

str(a)
matchLRPatterns("AAAAA", "CCCCC", 6, dna_string)


countPattern("ATTGCGCGA", dna_string, max.mismatch=2)
matchPattern("ATTGCGCGA",dna_string,max.mismatch=2)
