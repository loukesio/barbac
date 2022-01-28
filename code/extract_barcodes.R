################################
# libraries 
################################

library(here)         # set the working directory 
library(tidyverse)    # best package ever
library(Biostrings)   # read sequence data; be careful of the overlap with tidyverse
library(Rsamtools)

###############################
# set the working directory
###############################

setwd(here("data"))

###############################
# read data 
###############################
####

#$ the data contains 4 sequences that have been mapped to the reference

# extract the sequnce after using str_extract

readDNAStringSet("covert_bam_fasta.fasta") %>%   
  reverseComplement()   %>%     # I do not understand why we need the reverse complement 
  as.data.frame() %>%
  as_tibble() %>% 
  mutate(names=names(readDNAStringSet("covert_bam_fasta.fasta"))) %>% 
 dplyr::rename(seq=x) %>% 
  mutate(seq=str_extract(string = .$seq, pattern = "(?<=CTACG).*(?=CAGTC)")) %>% 
  mutate(length=str_length(seq))

library(tidystringdist)

df <- data.frame(stringsAsFactors = FALSE,
                 Name = as.factor(c(" CANON PVT. LTD ", " Antila,Thomas ", " Greg ",
                                    " St.Luke's Hospital ", " Z_SANDSTONE COOLING LTD ",
                                    " St.Luke's Hospital ", " CANON PVT. LTD. ",
                                    " SANDSTONE COOLING LTD ", " Greg ", " ANTILA,THOMAS ")),
                 City = as.factor(c(" Georgia ", " Georgia ", " Georgia ", " Georgia ",
                                    " Georgia ", " Georgia ", " Georgia ", " Georgia ",
                                    " Georgia ", " Georgia "))
)
df

df %>% 
  tidy_comb_all(Name) %>% 
  tidy_stringdist() %>% 
  filter(soundex == 0) %>% 
  pivot_longer(cols=starts_with("V"), 
               names_to = "match",
               values_to = "V"
  )

library("RBioinf")
Lpattern <- "CTACG"
Rpattern <- c("CAGTC")


x <- c("T", "A", "C", "G")

dict <- expand.grid(rep(list(x), 5)) %>% 
  unite("dictionary", 1:5, sep="")
dict

library(levenR)

L <- leven(dict$dictionary,Lpattern,substring2 = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::rename(distance="V1")  %>% 
  mutate(Lpattern=dict$dictionary) %>% 
  filter(distance==1 | distance==0) %>% 
  arrange(distance)
  
R <- leven(dict$dictionary,Rpattern,substring2 = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::rename(distance="V1")  %>% 
  mutate(Rpattern=dict$dictionary) %>% 
  filter(distance==1 | distance==0) %>%
  arrange(distance) %>% 
  select(Rpattern)

tidy_comb_all(Lpattern,Rpattern) 
  
    
  
setnames(.,refer) 
  mutate(barcodes= oper$barcodes) %>% 
  relocate(barcodes) %>% 

tidy_comb_all(dict, Lpattern)



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
