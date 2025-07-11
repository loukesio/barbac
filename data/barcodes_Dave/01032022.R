####################
# adding libraries
###################

library(tidyverse)    # best package for manipulating data
library(Biostrings)   # package for Bioinformatic analysis
library(stringdist)   # find the distance among strings
library(lvec)         # a cool package for saving memory
library(reclin2)      # Functions to assist in performing probabilistic record linkage and deduplication.
library(here)
library(data.table)
# add

##############################
# setting working directory
##############################

setwd(here("data","barcodes_Dave"))

#################
# read the data
#################

check <- c("CCACGCTGTCGTCCAACTCACGCAATCCATCGCGCAGTTTGTTCGGCAGGCTTTGTTAGTTC")
check %>% 
  str_length()

data <- fread(file = "4853_K barcode positions.txt")

sum(!complete.cases(data))


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
  dplyr::slice_sample(n=90000) %>% 
  mutate(barcode=as.character(barcode)) %>% 
  ungroup() %>% 
  mutate(id=row_number())

test
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

devtools::install_github("djvanderlaan/reclin2")

my_dist <- function(x, y, ...){
  1 - stringdist(x, y, method=c("lv"))
}

library(stringdist)
system.time(
p <- pair_minsim(test, on = "barcode", deduplication = TRUE,
                 default_comparator = my_dist, minsim = -2) # minsim=-2, reflects lv.distace = 3
)


 p2 <- p %>% 
  as.data.frame() %>% 
  mutate(across(1:2, .names="barcodes{col}", ~test$barcode[match(., test$id)])) %>% 
  mutate(across(4:5, .names="{col}.counts", ~ data$count[match(., data$barcode)]))

p2

data %>% 
  filter(grepl("ATTACGACTAGGTGAAAGCA", barcode))

head(data)

step4 <- inner_join(step3,step1) %>% 
  group_by(barcodes_ref) %>% 
  summarise(adds = sum(counts)) %>% 
  arrange(barcodes_ref)


add.counts <- inner_join(cluster1,oper) %>% 
  group_by(REFERENCE) %>% 
  summarise(adds = sum(n)) %>% 
  arrange(REFERENCE)
  


group_by(barcodes_ref) %>% 
  summarise(adds = sum(counts)) %>% 
  arrange(barcodes_ref)



######### Thats a big problem to solve ################

library(tidyverse)

df <- tibble(col1 = c("apple","apple","pple", "banana", "banana","bananna"),
             col2 = c("pple","app","app", "bananna", "banan", "banan"), 
             counts_col1 = c(100,100,2,200,200,2),
             counts_col2 = c(2,50,50,2,20,20),
             id=c(1,1,1,2,2,2))

df %>%
  rename_with(~ str_c("fruit_", .x), starts_with('col')) %>% 
  pivot_longer(cols = -id, names_to = c(".value", "grp"), 
               names_pattern = "(.*)_(col\\d+)") %>% 
  group_by(id, grp) %>%
  add_count(fruit) %>%
  group_by(id) %>% 
  summarise(central_fruit = fruit[which.max(n)], 
            fruits = toString(unique(fruit)), 
            sum_counts = sum(unique(counts)),
            counts = toString(sort(unique(counts), decreasing = TRUE)),
            .groups = 'drop' ) %>%
  relocate(counts, .before = 'sum_counts')


# add the missing step1 
# and the missing step3

# what is really missing from step3 
# It is missing to match hundred thousands of sequences to a reference 
library(reclin2)
library(stringdist)


reference = 

amatch(c("hello","g'day"),c("hi","hallo","ola"),maxDist=2)
amatch("leia",c("uhura","leela"),maxDist=5)

d <- data.frame(x=c(1,1,1,2,2,3,3), y=c(1,2,3,4,5,6,7), w=1:7)
d
# One-to-one matching:
d[match_n_to_m(d$x, d$y, d$w, n=999, m=999), ]


df1 <- data.frame(reference=c("cat","dog"))
df2 <- data.frame(data=c("cat","car","catt","cart","dog","dog","pitbull"), 
                  counts )
df1
df2
