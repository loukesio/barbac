# barcodes Dave 

setwd(here("data","barcodes_Dave"))

library(tictoc)

library(tidyverse)
library(Biostrings)
library(levenR)
library(data.table)
library(systemfonts)

##### set the plots style 

match_font('Avenir', italic = TRUE)

# test that your font has been indeed installed
ggplot(NULL, aes(0, 0)) +
  geom_text(
    aes(label = "The Avenir font"),
    size = 5, family = "Avenir"
  )

# change global theme setting for all following plots
theme_set(theme_bw(base_size = 20, base_family = "Avenir"))

# modify plot elements for all global plots
theme_update(
  plot.title = element_text(hjust = 0.5,size = 20, face = "bold", family="Avenir"), 
  legend.title.align=0.5,
  #panel.grid.major  = element_line(colour = "#F0F0F2", size=0.5),
  axis.ticks.x = element_line(colour = "#333333"),
  axis.ticks.y = element_line(colour = "#333333"),
  axis.ticks.length =  unit(0.26, "cm"), 
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  strip.background =element_rect(fill="white", colour="#333333")
)

####### here is not the plots style

tic("group barcodes")

data <- fread(file = "4853_K barcode positions.txt")

# how many 
dat1 <-  data %>% 
  group_by(corrected_insertion_site) %>% 
  count(corrected_insertion_site) %>% 
   arrange(desc(n)) 

dat1 %>% 
  ungroup() %>% 
  slice(-1) %>% 
  mutate(position= as.numeric(corrected_insertion_site)) %>% 
  ggplot(.) +
  geom_histogram(aes(log10(n), fill=n), binwidth = 1) +
  labs(x="barcode counts in log10 scale", y="counts", title="Distribution of barocde counts per position in log scale") +
  ggtext::geom_textbox(
  data = tibble(x = 3, y = 100000, label = "In the Plasmid Sequence, plasmidSeq, there are 41492 unique barcodes. The total number of barcodes in the plasmid is 129935."),
  aes(x, y, label = label),
  size = 4.5, family = "Avenir",
  fill = "cornsilk", box.color = "cornsilk3",
  width = unit(11, "lines")
)

# Barcode counts in the plasmidSeq
data %>% 
  filter(corrected_insertion_site == "plasmidSeq") %>% 
  summarise(sum(count))

test.distance <- data %>% 
  filter(corrected_insertion_site != "plasmidSeq" & corrected_insertion_site != "int64") %>% 
  dplyr::rename(position=corrected_insertion_site) %>% 
  group_by(position) %>% 
  filter(position == "1534")  

test.distance %>% 
  View()

# check for 1534 the following two barcodes the TTTTTATTTTACAGTTTTTT [83]
# and the TTTCCTCCAATAGACTTCTT [23]

levenAlign("TTTTTATTTTACAGTTTTTT","TTTCCTCCAATAGACTTCTT")

test.distance
leven(test.distance$barcode, substring2 = TRUE) %>% 
  as.data.frame() %>% 
  setnames(.,test.distance$barcode) %>% 
  mutate(barcodes=test.distance$barcode) %>%
  relocate(barcodes) %>% 
  pivot_longer(!barcodes, names_to = "barcodes_2", values_to = "distance") %>% 
  filter(distance>0 & distance <4)  %>% 
  unique(.[,c('barcodes','barcodes_2')])

  View()
  mutate(position=test.distance$position) 
  mutate(gene_name=test.distance$gene_name) %>% 
  mutate(count=test.distance$count) %>% 
  mutate(total = sum(count)) %>% 
  filter(count == max(count)) 




test <- data %>% 
  filter(corrected_insertion_site != "plasmidSeq" & corrected_insertion_site != "int64") %>% 
dplyr::rename(position=corrected_insertion_site) %>% 
  filter(position == "1534") %>% 
  group_by(position) %>% 
  do(lv.dist= {
    check <- leven(test$barcode, substring2 = TRUE) %>% 
      as.data.frame() %>% 
      setnames(.,test$barcode) %>% 
      mutate(barcodes=test$barcode) %>%
      relocate(barcodes) %>% 
      mutate(position=test$position) %>% 
      mutate(gene_name=test$gene_name) %>% 
      mutate(count=test$count) 
    check
  })





# install.packages("tabulate")
# library(tabulate)
# department <- tabulate_table() %>% 
#   table_add_row(c("Research", "Engineering"))
# department

test
test$lv.dist

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


