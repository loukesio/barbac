################################
# libraries 
################################

library(here)         # set the working directory 
library(tidyverse)    # best package ever
library(Biostrings)   # read sequence data; be careful of the overlap with tidyverse
library(levenR)       # very useful package 
library(tidystringdist)

###############################
# set the working directory
###############################

setwd(here("data"))

# First I am creating a dictionary  of 5 base pair bases

x <- c("T", "A", "C", "G")
dict <- expand.grid(rep(list(x), 5)) %>% 
  unite("dictionary", 1:5, sep="")

Lpattern <- c("CTACG")
Rpattern <- c("CAGTC")

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

crossing(L,R) # is a valid option and then you can achieve what you want

#___________________________________________________________________________________________________________
# ╦═╗╔═╗╔═╗╔═╗╦═╗╔═╗╔╗╔╔═╗╔═╗╔═╗
# ╠╦╝║╣ ╠╣ ║╣ ╠╦╝║╣ ║║║║  ║╣ ╚═╗
# ╩╚═╚═╝╚  ╚═╝╩╚═╚═╝╝╚╝╚═╝╚═╝╚═╝

# 1. https://stackoverflow.com/questions/11388359/unique-combination-of-all-elements-from-two-or-more-vectors [how to make all combinations]

