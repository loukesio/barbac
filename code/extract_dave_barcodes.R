# barcodes Dave 

library(tictoc)
library(here)
library(tidyverse)
library(Biostrings)
library(levenR)
library(data.table)
library(systemfonts)

##### set the plots style 

setwd(here("data","barcodes_Dave"))

########################################

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

#___________________________________________
# ┌┬┐┌─┐┌┬┐┌─┐┬    ┌┐ ┌─┐┬─┐┌─┐┌─┐┌┬┐┌─┐┌─┐
#  │ │ │ │ ├─┤│    ├┴┐├─┤├┬┘│  │ │ ││├┤ └─┐
#  ┴ └─┘ ┴ ┴ ┴┴─┘  └─┘┴ ┴┴└─└─┘└─┘─┴┘└─┘└─┘
#___________________________________________

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
  labs(x="barcode counts in log10 scale", y="counts", title="Distribution of barocde counts in log scale") +
  ggtext::geom_textbox(
  data = tibble(x = 2.5, y = 120000, label = "In the Plasmid Sequence, plasmidSeq, there are 41492 unique barcodes. The total number of barcodes in the plasmid is 129935."),
  aes(x, y, label = label),
  size = 4.5, family = "Avenir",
  fill = NA, box.color = NA,
  width = unit(11, "lines")
) +
  geom_point(data=tibble(x1=2.5, y1=125000), aes(x1,y1), size=90, shape=21)

#______________________________________________________
# ┌┐ ┌─┐┬─┐┌─┐┌─┐┌┬┐┌─┐┌─┐  ┬┌┐┌  ┌─┐┬  ┌─┐┌─┐┌┬┐┬┌┬┐
# ├┴┐├─┤├┬┘│  │ │ ││├┤ └─┐  ││││  ├─┘│  ├─┤└─┐││││ ││
# └─┘┴ ┴┴└─└─┘└─┘─┴┘└─┘└─┘  ┴┘└┘  ┴  ┴─┘┴ ┴└─┘┴ ┴┴─┴┘
#______________________________________________________

data %>% 
  filter(corrected_insertion_site == "plasmidSeq") %>% 
  summarise(sum_plasmid=sum(count))

#_________________________________________________________________________________
# ╔╗ ┌─┐┬─┐┌─┐┌─┐┌┬┐┌─┐  ╦╔╗╔  ╔╦╗┬┌─┐┌─┐┌─┐┬─┐┌─┐┌┐┌┌┬┐  ╔═╗┌─┐┌─┐┬┌┬┐┬┌─┐┌┐┌┌─┐
# ╠╩╗├─┤├┬┘│  │ │ ││├┤   ║║║║   ║║│├┤ ├┤ ├┤ ├┬┘├┤ │││ │   ╠═╝│ │└─┐│ │ ││ ││││└─┐
# ╚═╝┴ ┴┴└─└─┘└─┘─┴┘└─┘  ╩╝╚╝  ═╩╝┴└  └  └─┘┴└─└─┘┘└┘ ┴   ╩  └─┘└─┘┴ ┴ ┴└─┘┘└┘└─┘
#_________________________________________________________________________________

head(data)
diff.positions <- data %>% 
  filter(corrected_insertion_site != "plasmidSeq") %>% 
  dplyr::rename(position=corrected_insertion_site) %>% 
  group_by(barcode) %>% 
  summarise(groups = list(position), n_groups = n(), counts = sum(count)) %>% 
  arrange(desc(n_groups)) %>% 
  filter(n_groups > 1) %>% 
  select(-groups) %>% 
  mutate(n_positions=n_groups)

data %>% 
  filter(corrected_insertion_site != "plasmidSeq") %>% 
  dplyr::rename(position=corrected_insertion_site) %>% 
  group_by(barcode) %>% 
  summarise(groups = list(position), n_groups = n(), counts = sum(count)) %>% 
  arrange(desc(n_groups)) %>% 
  filter(n_groups < 2) %>% 
  filter(grepl("^GGG", barcode))

df1 %>% filter(!grepl("^x|xx$", fruit))

  grepl('^GG', x) # starts with AB?


write.table(diff.positions,"barcodes_diffpos.txt")

  ungroup() %>% 
  count(barcode)
  # slice(1) %>% 
  # select(groups) %>% 
  # unlist()

#_________________________________________________________
# ╔═╗┬  ┬ ┬┌─┐┌┬┐┌─┐┬─┐  ┌─┐┌─┐┬─┐  ╔═╗┌─┐┌─┐┬┌┬┐┬┌─┐┌┐┌
# ║  │  │ │└─┐ │ ├┤ ├┬┘  ├─┘├┤ ├┬┘  ╠═╝│ │└─┐│ │ ││ ││││
# ╚═╝┴─┘└─┘└─┘ ┴ └─┘┴└─  ┴  └─┘┴└─  ╩  └─┘└─┘┴ ┴ ┴└─┘┘└┘
#_________________________________________________________

data %>% 
  filter(corrected_insertion_site != "plasmidSeq" & corrected_insertion_site != "int64") %>% 
  dplyr::rename(position=corrected_insertion_site) %>% 
  group_by(position) %>% 
  filter(position == "1534")  %>% 
  arrange(desc(count))
  

# check for 1534 the following two barcodes the TTTTTATTTTACAGTTTTTT [83]
# and the TTTCCTCCAATAGACTTCTT [23]

levenAlign("TTTTTATTTTACAGTTTTTT","TTTCCTCCAATAGACTTCTT")

test.distance
leven(test.distance$barcode, substring2 = TRUE) %>% 
  as.data.frame() %>% 
  setnames(.,test.distance$barcode) %>% 
  mutate(barcode=test.distance$barcode) %>%
  relocate(barcode) %>% 
  pivot_longer(!barcode, names_to = "barcodes_2", values_to = "distance") %>% 
  filter(distance>0 & distance <4) %>% 
    left_join(.,test.distance, by="barcode") %>%
  mutate(Var = map2_chr(barcode, barcodes_2, ~toString(sort(c(.x, .y))))) %>%
  distinct(Var, .keep_all = TRUE) %>%
  select(-Var) %>%
  arrange(barcode) %>% 
  left_join(., test.distance, by = c("barcodes_2" = "barcode")) %>% 
  select(barcode,barcodes_2, distance, position.x, gene_name.x, count.x, count.y) %>% 
  rename(count_barcode=count.x,count_barcodes_2=count.y)
  

test.distance

  unique(.[,c('barcodes','barcodes_2')])

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


