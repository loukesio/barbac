################################
# libraries 
################################

library(here)            # set the working directory 
library(tidyverse)       # best package ever
library(Biostrings)      # read sequence data; be careful of the overlap with tidyverse
library(levenR)          # very useful package 
library(tidystringdist)  # finding string distance an alternative to levenR

###############################
# set the working directory
###############################

setwd(here("data"))

# First I am creating a dictionary  of 5 base pair bases

#____________________________________________________
# ┌─┐┌─┐┌┬┐┌┬┐┌─┐┬─┐┌┐┌┌─┐  ┌┬┐┌─┐┌┬┐┌─┐┌┐ ┌─┐┌─┐┌─┐
# ├─┘├─┤ │  │ ├┤ ├┬┘│││└─┐   ││├─┤ │ ├─┤├┴┐├─┤└─┐├┤ 
# ┴  ┴ ┴ ┴  ┴ └─┘┴└─┘└┘└─┘  ─┴┘┴ ┴ ┴ ┴ ┴└─┘┴ ┴└─┘└─┘
#____________________________________________________

x <- c("T", "A", "C", "G")
dict <- expand.grid(rep(list(x), 5)) %>% 
  unite("dictionary", 1:5, sep="")

Lpattern <- c("CTACG")
Rpattern <- c("CAGTC")

L <- leven(dict$dictionary,Lpattern,substring2 = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::rename(distL="V1")  %>% 
  mutate(Lpattern=dict$dictionary) %>% 
  filter(distL ==1 | distL==0) %>% 
  arrange(distL)

R <- leven(dict$dictionary,Rpattern,substring2 = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::rename(distR="V1")  %>% 
  mutate(Rpattern=dict$dictionary) %>% 
  filter(distR==1 | distR==0) %>%
  arrange(distR) 

index.dict <- crossing(L,R) # is a valid option and then you can achieve what you want

# pattern1 database
patterns1 <- paste0(index.dict$Lpattern, "(.*)", index.dict$Rpattern)
names(patterns1) <- paste0("pattern", 1:length(patterns1))

# pattern2 database
patterns2 <- paste0("(?<=", index.dict$Lpattern, ")(.*)(?=", index.dict$Rpattern, ")")
names(patterns2) <- paste0("pattern", 1:length(patterns2))

#___________________________________________________
# ╔═╗═╗ ╦╔╦╗╦═╗╔═╗╔═╗╔╦╗  ╔╗ ╔═╗╦═╗╔═╗╔═╗╔╦╗╔═╗╔═╗  
# ║╣ ╔╩╦╝ ║ ╠╦╝╠═╣║   ║   ╠╩╗╠═╣╠╦╝║  ║ ║ ║║║╣ ╚═╗  
# ╚═╝╩ ╚═ ╩ ╩╚═╩ ╩╚═╝ ╩   ╚═╝╩ ╩╩╚═╚═╝╚═╝═╩╝╚═╝╚═╝  
#___________________________________________________

# first read the sequence from which you want to extract barcodes 

raw.data <- readDNAStringSet("covert_bam_fasta.fasta") %>%   
  reverseComplement() %>% 
  as.data.frame() %>%
  as_tibble() %>% 
  mutate(names=names(readDNAStringSet("covert_bam_fasta.fasta"))) %>% 
  dplyr::rename(seq=x)

test <- expand_grid(seq = raw.data$seq,
            pattern = patterns2) %>%
  distinct() %>% 
  mutate(match = str_extract_all(seq, pattern)) %>% 
  pivot_wider(
    names_from = "pattern",
    values_from = "match"
  ) %>% 
  dplyr::rename_with(~names(patterns2),
              .cols = -seq) %>% 
  select(where(~!all(lengths(.) == 0))) %>% 
  mutate(across(starts_with('pattern'), ~ lapply(., function(x) ifelse(length(x) == 0, NA, x))))

test

# i ll find a way   
test1 <- test %>%   
as.data.frame() %>% 
mutate_at(vars(2:last_col()), as.character) %>% 
pivot_longer(!seq,names_to = "patterns", values_to = "barcodes") %>% 
  arrange(patterns) %>% 
  mutate(barcode_length=str_length(barcodes)) %>% 
  filter(barcode_length < 27 & barcode_length > 23)

test1

test1 %>% 
  group_by(barcodes) %>% 
  mutate(count=n())
# how to rename columns which have lists inside and then how to select only the  rows from the ones that are not NA


#___________________________________________________________________________________________________________
# ╦═╗╔═╗╔═╗╔═╗╦═╗╔═╗╔╗╔╔═╗╔═╗╔═╗
# ╠╦╝║╣ ╠╣ ║╣ ╠╦╝║╣ ║║║║  ║╣ ╚═╗
# ╩╚═╚═╝╚  ╚═╝╩╚═╚═╝╝╚╝╚═╝╚═╝╚═╝

# 1. https://stackoverflow.com/questions/11388359/unique-combination-of-all-elements-from-two-or-more-vectors [how to make all combinations]

library(ggplot2)
library(paint)
library(here)
library(data.table)
library(tidyverse)
library(Biostrings)
library(rBLAST)
library(reprex)
library(styler)
setwd(here("data","test_extraction"))

###############
# Reference 
#############

readDNAStringSet("Reference_barcodes.fasta") %>% 
  as.data.frame()

library(GenomicAlignments)

readDNAStringSet("")

stack1 <- stackStringsFromBam("un_mapped.sorted.bam", 
                             param=GRanges("Reference_barcodes:54-78"))
stack1


stack2 <- stackStringsFromBam("unmapped.alignment.sorted.bam", 
                              param=GRanges("Reference_barcodes:54-78"))


# https://wikis.utexas.edu/display/bioiteam/Adva
# use this for tomorrow as well https://stackoverflow.com/questions/4854130/how-to-iterate-over-file-names-in-a-r-script

stack
  (alphabetFrequency(., as.prob = T,baseOnly=T))

afmc=consensusMatrix(stack, baseOnly=T,as.prob = T)
tafmc=t(afmc)
matplot(tafmc[,-5], type="l", lwd=2, xlab="Read Length", ylab= "Base frequency at each position")
legend(legend = colnames(tafmc)[-5],"topright",col=1:4, lty=1:4, lwd=2)


# essential adding in the plot
