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



library(stringdist)
stringdist("ca","abc")

stringdist(c("ATTTC","TTTC"), c("ATTTC","TTTC"),method="lv", nthread = 4)

a=c("ATTTC","TTTC")

data(phiX174Phage)
phiX174Phage


stringDist(phiX174Phage, method = "levenshtein", nthread=2)

library(stringdist)
vec = c("keep", "teem", "meat", "weep")

out <- as.matrix(stringdistmatrix(vec), method=levenshtein)
dimnames(out) <- list(vec, vec)
out

# References 

# #https://jamesfeigenbaum.github.io/post/2017/speed-testing-string-distance-measures/
library(tidyverse)

df <- data.frame(names=c("A ADAM", "S BEAN", "A APPLE", "J BOND", "J BOND"), 
                 v1=c("Test_a", "Test_b", "Test_a", "Test_b", "Test_b"), 
                 v2=c("Test_c", "Test_c", "Test_d", "Test_d", "Test_d"))

df

lapply(df, stringdist::stringdistmatrix, method = 'lv') %>% 
   map_df(broom::tidy, .id = 'var') %>% 
  spread(var, distance)


library(rbenchmark)
library(RecordLinkage)
library(stringdist)

 x <- sapply(sample(5:25,1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse="")) 
 y <- sapply(sample(5:25,1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse="")) 

 benchmark(
     d1 <- stringdist(x,y,method='lv'),
      d2 <- levenshteinDist(x,y),
     replications=10)
