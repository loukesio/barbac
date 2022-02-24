# barcodes Dave 

setwd(here("data","barcodes_Dave"))

library(tictoc)

library(tidyverse)
library(Biostrings)
library(levenR)

tic("clustering")

step1 <- readDNAStringSet("example.fasta") %>% 
  as.data.frame() %>%
  as_tibble() %>% 
  dplyr::slice(-1) %>% 
  dplyr::rename(barcodes=`x`) %>%
  dplyr::count(barcodes) %>%
  dplyr::arrange(desc(n)) %>% 
  dplyr::rename(counts=n)

step2 <- step1 %>% 
  dplyr::mutate(length=str_length(barcodes)) %>% 
  dplyr::filter(length > 19 & length  < 21) %>% 
  slice(1:10) %>%              ####### !!!!!!!!! remove it if necessary
  dplyr::select(barcodes) 

step3 <- leven(step2$barcodes,substring2 = TRUE)  %>% 
  as.data.frame() %>% 
  set_names(.,step2$barcodes) %>% 
  mutate(barcodes_ref = step2$barcodes) %>% 
  relocate(barcodes_ref) %>% 
  pivot_longer(!barcodes_ref, names_to = "barcodes", values_to = "dist") %>% 
  mutate(clust = as.numeric(forcats::as_factor(barcodes_ref))) %>%   # could also use data.table::rleid
  filter(dist<3 & dist>0)  %>% 
  arrange(clust)

step4 <- inner_join(step3,step1) %>% 
  group_by(barcodes_ref) %>% 
  summarise(adds = sum(counts)) %>% 
  arrange(barcodes_ref)

#write.table(step4,"leven_example.txt", sep=",")

toc()


### alternative clustering methods

# here are the things to do 


# I think you have already work for this 
# you find the counts for both barcodes # you can do that with the join, I am sure you have worked on that on the step2
# and then you have found the function to keep the highest 
# https://stackoverflow.com/questions/71161281/sum-the-counts-of-the-unique-elements-of-a-group-and-keep-the-one-with-the-highe
# but also break the 

