# barcodes Dave 

library(here)
library(tidyverse)
library(Biostrings)
library(levenR)
library(data.table)    # for using fread
library(systemfonts)   # add teh fonts that you like
library(stringr)       # this is important for reading the length of the strings in R
library(patchwork)     # combine ggplots together
library(ltc)           # library for coloring projects
library(ggtext)        # simple Rmarkdown and HTML rendering for ggplot2
library(vwr)           # install vwr library for string clustering with Levenstein distance
library(xlsx)          # read excel files
library(viridis)
library(vegan)
library(betapart)
library(ggrepel)


setwd(here("data","barcodes_Dave"))


step1 <- readDNAStringSet("example.fasta") %>% 
  as.data.frame() %>%
  as_tibble() %>% 
  dplyr::slice(-1) %>% 
  dplyr::rename(barcodes=`x`) %>%
  dplyr::count(barcodes) %>%
  dplyr::arrange(desc(n)) %>% 
  rename(counts=n)

step2 <- step1 %>% 
  mutate(length=str_length(barcodes)) %>% 
  filter(length > 19 & length  < 21) %>% 
  slice(1:10)  %>% 
  select(barcodes)
  
step3 <- leven(step2$barcodes,substring2 = TRUE)  %>% 
  as.data.frame() %>% 
  set_names(.,step2$barcodes) %>% 
  mutate(barcodes2 = step2$barcodes) %>% 
  relocate(barcodes2) %>% 
  pivot_longer(!barcodes2, names_to = "barcodes1", values_to = "distance")

step3
  filter_at(vars(2:11), any_vars(. < 10 & .>0)) 
step3


cluster1 <- leven(oper$barcodes,refer,substring2 = TRUE) %>% 
  as.data.frame() %>% 
  setnames(.,refer) %>% 
  mutate(barcodes= oper$barcodes) %>% 
  relocate(barcodes) %>% 
  pivot_longer(-barcodes, names_to = "REFERENCE", values_to = "dist") %>% 
  mutate(clust = as.numeric(forcats::as_factor(REFERENCE))) %>%  # could also use data.table::rleid
  filter(dist<3 & dist>0)  %>% 
  relocate(REFERENCE) %>% 
  arrange(desc(clust)) 

step1 %>% 
  leven(.$barcodes,substring2 = TRUE)



  stringdistmatrix(barcodes, useNames=TRUE ,method = "lv") 
  as.matrix() %>% 
  heatmap()




