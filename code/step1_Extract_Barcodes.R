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
            pattern = patterns1) %>%
  distinct() %>% 
  mutate(match = str_extract_all(seq, pattern)) %>% 
  pivot_wider(
    names_from = "pattern",
    values_from = "match"
  ) %>% 
  dplyr::rename_with(~names(patterns1),
              .cols = -seq) %>% 
  select(where(~!all(lengths(.) == 0))) %>% 
  as.data.frame() 
  #mutate_all(., function(x) lapply(x, length) == 0) <- NA )

test  

test[apply(test, 2, function(x) lapply(x, length) == 0)] <- NA #Replace empty lists with NA

test %>% 
  as.data.frame() %>% 
  mutate_at(vars(2:last_col()), as.character) %>% 
  pivot_longer(!seq,names_to = "patterns", values_to = "barcodes")




test[apply(test, 2, function(x) lapply(x, length) == 0)] <- NA #Replace empty lists with NA

test %>% 
  as.data.frame()

check = (test[apply(test, 2, function(x) lapply(x, length) == 0) == T] <- NA) #Replace empty lists with NA
check
test %>% 
  unnest(vars(2:last_col()), keep_empty = TRUE)

test1
test
test %>% 
  mutate_at(vars(2:last_col()), funs(unlist))
  
  

select(pattern1:pattern2) %>% 
  unnest(pattern1)
# how to rename columns which have lists inside and then how to select only the  rows from the ones that are not NA


#___________________________________________________________________________________________________________
# ╦═╗╔═╗╔═╗╔═╗╦═╗╔═╗╔╗╔╔═╗╔═╗╔═╗
# ╠╦╝║╣ ╠╣ ║╣ ╠╦╝║╣ ║║║║  ║╣ ╚═╗
# ╩╚═╚═╝╚  ╚═╝╩╚═╚═╝╝╚╝╚═╝╚═╝╚═╝

# 1. https://stackoverflow.com/questions/11388359/unique-combination-of-all-elements-from-two-or-more-vectors [how to make all combinations]

