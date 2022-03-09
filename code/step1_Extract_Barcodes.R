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
setwd(here("data"))


library(Rsamtools)
bam <- scanBam("test_4_seqs_sorted.bam")
bam

#names of the BAM fields
names(bam[[1]])

#distribution of BAM flags
table(bam[[1]]$flag)

.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
bam_df
names(bam_df) <- bam_field

bam_df 


bam_df %>% 
  as_tibble() %>% 
  select(seq,cigar) %>% 
  View()
bam_df
test <- bam_df %>% 
  as_tibble() %>% 
  select(seq,flag,cigar) %>% 
separate(cigar, into = str_c("col", 1:5), 
         sep = "(?<=\\D)(?=\\d)", fill = "right", remove = FALSE) %>% 
mutate(
    across(starts_with("col"),~case_when(
      is.na(.) ~ NA_real_,
      grepl("[SDI]$", .) ~ parse_number(.),
      TRUE ~ 0
    )
    ))

test %>% 
  View()


test %>% 
  mutate(removed = str_sub(seq, col1 + 1)) %>% 
  select(removed) %>% 
  mutate(check=str_sub(removed, start=54, end=78)) %>% 
  View()

TACGGTTATATTGACAGACCGAGGG

