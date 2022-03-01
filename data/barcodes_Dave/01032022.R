####################
# adding libraries
###################

library(tictoc) 
library(tidyverse)
library(Biostrings)
library(stringdist)

##############################
# setting working directory
##############################

setwd(here("data","barcodes_Dave"))

#################
# read the data
################

data <- fread(file = "4853_K barcode positions.txt")

step1 <- data %>% 
  filter(corrected_insertion_site != "plasmidSeq") %>% 
  dplyr::rename(position=corrected_insertion_site) %>% 
  group_by(barcode) %>% 
  summarise(groups = list(position), n_groups = n(), counts = sum(count)) %>% 
  arrange(desc(n_groups)) %>% 
  filter(n_groups < 2) %>% 
  dplyr::mutate(length=str_length(barcode)) %>% 
  dplyr::filter(length > 19 & length  < 21) %>% 
  dplyr::select(barcode)

set.seed(123)
test <- step1 %>% 
  dplyr::slice_sample(n=10) %>% 
  mutate(barcode=as.character(barcode))


#_________________________
# ┌┬┐┌─┐┌┬┐┬ ┬┌─┐┌┬┐  ╔═╗
# │││├┤  │ ├─┤│ │ ││  ╠═╣
# ┴ ┴└─┘ ┴ ┴ ┴└─┘─┴┘  ╩ ╩
#_________________________

system.time(
  method.A <- stringdistmatrix(test$barcode, method=c("lv"), nthread = 16) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  set_names(.,test$barcode) %>% 
  mutate(barcodes_ref = test$barcode) %>% 
  relocate(barcodes_ref) %>% 
  pivot_longer(!barcodes_ref, names_to = "barcodes", values_to = "dist") %>% 
  mutate(clust = as.numeric(forcats::as_factor(barcodes_ref))) %>%   # could also use data.table::rleid
  filter(dist< 4 & dist>0)  %>% 
  arrange(clust)
)

#________________________
# ┌┬┐┌─┐┌┬┐┬ ┬┌─┐┌┬┐  ╔╗ 
# │││├┤  │ ├─┤│ │ ││  ╠╩╗
# ┴ ┴└─┘ ┴ ┴ ┴└─┘─┴┘  ╚═╝
#________________________


my_dist <- function(x, y, ...){
  1 - stringdist(x, y, method=c("lv"))
}

system.time(
p <- pair_minsim(test, on = "barcode", deduplication = TRUE,
                 default_comparator = my_dist, minsim = -2)
)
p
