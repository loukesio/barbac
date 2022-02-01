library(here)           # add the directory
library(levenR)         # see the distance 
library(tidyverse)      # best package ever
library(Biostrings)     # analyse fasta|q data
library(systemfonts)    # put a new font in your data
library(tidystringdist) # add distacen of strings in R

setwd(here("data"))
getwd()

# add a dictionary 

x <- c("T", "A", "C", "G")
dict <- expand.grid(rep(list(x), 5)) %>% 
  unite("dictionary", 1:5, sep="")
dict

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

database <- tidy_comb_all(L[,2],R) %>% 
  dplyr::rename(Lpattern=V1,Rpattern=V2)

database

database2=data.frame(L=c("CTAGG","CTCCG"), R=c("CAGTC","CAGTC"))
database2

CCACGAAGCTCTCCTACGTACGGTTATATTGACAGACCGAGGGCAGTCCAGCGCCAACCAGATAAGTGAAATCTAGTTCCA
CCACGAAGCTCTCCTACGTACGGTTATATTGACAGACCGAGGGCAGTCCAGCGCCAACCAGATAAGTGAAATCTAGTTCCA
CCACGAAGCTCTCCTAGGGGGGGGCTATTTTGGACTGCGTTACCAGTCCAGCGCCAACCAGATAAGTGGAATCTAGTTCGA
CCACGTAGCTCTCCTCCGTGCGGTTATATTGACAGACCGAGGGCAGTCCAGCGCCAACCAGATAAGTGAAATCTAGTTCCA



# here with strict border I am missing two from the output

readDNAStringSet("convert_test_special.fasta") %>%   
    reverseComplement() %>%        # I do not understand why we need the reverse complement 
    as.data.frame() %>%
    as_tibble() %>% 
    dplyr::rename(seq=x) %>% 
    mutate(barcodes=str_extract(string = .$seq, pattern = "(?<=CTACG).*(?=CAGTC)")) %>% 
    mutate(name=names(readDNAStringSet("convert_test_special.fasta"))) %>% 
  filter(is.na(barcodes)) %>% 
  mutate(barcodes2= str_extract(string=.$seq, 
  str_c("(?<=",str_c(database2[1,1], collapse = ""),").*(?=",str_c(database2[1,2], collapse = ""),")"))) %>% 
  View()
  

patterns <- paste0(database$L, "(.*)", database$R)

names(patterns) <- paste0("pattern", 1:3)

cbind(
  data,
  lapply(
    patterns,
    \(x) str_match(data$seq, x)[,2]
  )
)












database[1,1]
str_c(database[1,1], collapse = "")

pattern2 <- str_c("(?<=",str_c(database[2,1], collapse = ""),").*(?=",str_c(database[2,2], collapse = ""),")")
pattern2


check1 <- readDNAStringSet("convert_test_special.fasta") %>%   
  reverseComplement() %>%        # I do not understand why we need the reverse complement 
  as.data.frame() %>%
  as_tibble() %>% 
  dplyr::rename(seq=x) %>% 
  mutate(name=names(readDNAStringSet("convert_test_special.fasta"))) %>% 
  mutate(barcodes= str_extract(string = .$seq, pattern = "(?<=CTACG).*(?=CAGTC)"))

check1 %>% 
  filter_all(any_vars(is.na(.))) %>% 
  mutate(barcodes=str_extract())





check2 <- check1 %>% 
  filter(is.na(barcodes)) 

require("dplyr")
Sequences=c("AzzY","BbDe")
Database=c("TTUAzzY","aaa","DBbDe","CAzzY")

df=as.data.frame(sapply(Sequences, function(x) grepl(x,Database)))
Sequences
stats=df %>% summarise_each(funs(sum))
stats
cbind(Sequences,as.numeric(stats))




delim1 <- "TGCA"
delim2 <- "CTGC"

for (i in 0:nchar(delim1))
{
  for (j in 0:nchar(delim2))
  {
    sdelim1 <- str_split(delim1, "") %>% unlist
    sdelim2 <- str_split(delim2, "") %>% unlist
    
    if (i != 0)
      sdelim1[i] <- "[A-Z]"
    
    if (j != 0)
      sdelim2[i] <- "[A-Z]"
    
    pattern <- str_c("(?<=",str_c(sdelim1, collapse = ""),").*(?=",str_c(sdelim2, collapse = ""),")")
    
    print(str_extract(s,pattern))
    
  }
}
