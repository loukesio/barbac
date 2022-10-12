#___________________________________________
#  _  _  _                     _          
# | |(_)| |__  _ _  __ _  _ _ (_) ___  ___
# | || || '_ \| '_|/ _` || '_|| |/ -_)(_-<
# |_||_||_.__/|_|  \__,_||_|  |_|\___|/__/
#                                         
#___________________________________________

# here are the libraries that I need 
# what is the meaning of the lib.loc? 
# lib.loc is a character vector describing the location of R library trees to search through, or NULL 
library(tidyverse, lib.loc="/data/modules/R/4.1.2/lib64/R/library")
library(data.table)                                                    # best package together with tidyverse for data analysis
library(here)
library(gt)
library(gtExtras)
library(stringdist)                                                    # best package for calculating string distances
library(reclin2, lib.loc="/data/modules/R/4.1.2/lib64/R/library")      # package to minimize memory in R

#______________________________________________
# ┬ ┬┌─┐┬─┐┬┌─┬┌┐┌┌─┐  ┌┬┐┬┬─┐┌─┐┌┬┐┌─┐┬─┐┬ ┬
# ││││ │├┬┘├┴┐│││││ ┬   │││├┬┘├┤  │ │ │├┬┘└┬┘
# └┴┘└─┘┴└─┴ ┴┴┘└┘└─┘  ─┴┘┴┴└─└─┘ ┴ └─┘┴└─ ┴ 
#______________________________________________

setwd("/home/theodosiou/Projects/Barcodes/Static_Experiment/Amplicon_Barcodes/rep2/merged/assembled_reads/minimap2_output/") # setting the working directory, exactly where you want it

#___________________________________
# ╔═╗╔═╗╦═╗  ╔═╗╔╗╔╔═╗  ╔═╗╦╦  ╔═╗
# ╠╣ ║ ║╠╦╝  ║ ║║║║║╣   ╠╣ ║║  ║╣ 
# ╚  ╚═╝╩╚═  ╚═╝╝╚╝╚═╝  ╚  ╩╩═╝╚═╝
#__________________________________

files=list.files(pattern = "*barcodes.txt$")   # select all the files that have isolate barcodes 

day1= files[1]                                 # From all these files work only with the one from the first day --> aka "day1"

#___________________________________
# ┌─┐┬ ┬┌┐┌┌─┐┌┬┐┬┌─┐┌┐┌  ╔═╗╔╗╔╔═╗
# ├┤ │ │││││   │ ││ ││││  ║ ║║║║║╣             # read raw file and export 📊 statistic out it
# └  └─┘┘└┘└─┘ ┴ ┴└─┘┘└┘  ╚═╝╝╚╝╚═╝
#___________________________________

files=list.files(pattern = "*fastq.bam$") 
files

# for(i in 1:length(files)) {
#   
#   stack <- stackStringsFromBam(files[i], 
#                                param=GRanges("Reference_barcodes:54-78"),
#                                use.names = TRUE,
#                                what = "seq") 
#   data <- stack %>% 
#     as.data.frame() %>% 
#     dplyr::rename(barcode=x) %>% 
#     group_by(barcode) %>% 
#     dplyr::count() %>% 
#     arrange(desc(n)) %>% 
#     dplyr::rename(counts=n) %>% 
#     mutate(length_barcode=str_length(barcode)) %>% 
#     filter(length_barcode < 27 | length_barcode > 23)
#   
#   write.table(data, quote=FALSE, sep=",", sub(".assembled.fastq_sorted.bam","_barcodes.txt", files[i]))
#   
# }



barcodes.extract.m1 <- function(file) {
  
  raw.barcodes <- read.table(file,sep=",") %>% 
    dplyr::mutate(barcode=str_replace_all(barcode,"-","")) %>% 
    dplyr::mutate(length_barcode=str_length(barcode)) %>% 
    filter(length_barcode!="0")
  
  raw.barcodes %>% 
    mutate(.,bin = case_when((length_barcode < 22 & length_barcode > 0)  ~ "(0,22)",
                             (length_barcode >=22 & length_barcode <= 28) ~ "[22,28]",
                             TRUE ~ "[29,∞)")) %>% 
    group_by(bin) %>% 
    count() %>% 
    as.data.frame() %>% 
    gt() %>% 
    tab_header(
      title = "day1"
    ) %>% 
    cols_label(
      n = "No of barcodes",
    ) %>% 
    tab_source_note(
      source_note = "Barbac: A versatile tool to quantify barcodes."
    ) %>% 
    gt_theme_espn() %>% 
    gt::gtsave(filename = here(paste0(str_replace(files[i], pattern = "S(\\d+_\\d+).+", replacement = "\\1_summary_rawdata"),".png")))
  
plot.rawbar <- ggplot(raw.barcodes, aes(x=length_barcode)) +
    geom_histogram(fill="#FF5349", color=NA, alpha=0.8) +
    scale_x_continuous(breaks=seq(0,30, by=5), 
                       limits=c(min(raw.barcodes$length_barcode)-5,max(raw.barcodes$length_barcode) +5)) +
    labs(x="Barcode length (bp)", y="Number of barcodes", 
         title=paste0(str_replace(files[i], pattern = "S(\\d+_\\d+).+", replacement = "\\1_raw_data")), caption=c("Barbac: A versatile tool to quantify barcodes")) + # find a way to make it more general
    theme_bw() 
  
  ggsave(here(paste0(str_replace(files[i], pattern = "S(\\d+_\\d+).+", replacement = "\\1_barcode_length"),".png")))
  
  
  return(raw.barcodes)
}

raw.day1 <- barcodes.extract.m1(day1)

#___________________________________
# ┌─┐┬ ┬┌┐┌┌─┐┌┬┐┬┌─┐┌┐┌  ╔╦╗╦ ╦╔═╗
# ├┤ │ │││││   │ ││ ││││   ║ ║║║║ ║
# └  └─┘┘└┘└─┘ ┴ ┴└─┘┘└┘   ╩ ╚╩╝╚═╝
#___________________________________
library(levenR)

barcodes.extract.m2 <- function(distance, file) {

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

index.dict <- crossing(L,R) # is a
head(index.dict)

patterns1 <- paste0(index.dict$Lpattern, "(.*)", index.dict$Rpattern)

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

return(test)

}


#___________________________________________________________________________________________
#
#           $$\                       $$\                         $$\                     
#           $$ |                      $$ |                        \__|                    
#  $$$$$$$\ $$ |$$\   $$\  $$$$$$$\ $$$$$$\    $$$$$$\   $$$$$$\  $$\ $$$$$$$\   $$$$$$\  
# $$  _____|$$ |$$ |  $$ |$$  _____|\_$$  _|  $$  __$$\ $$  __$$\ $$ |$$  __$$\ $$  __$$\ 
# $$ /      $$ |$$ |  $$ |\$$$$$$\    $$ |    $$$$$$$$ |$$ |  \__|$$ |$$ |  $$ |$$ /  $$ |
# $$ |      $$ |$$ |  $$ | \____$$\   $$ |$$\ $$   ____|$$ |      $$ |$$ |  $$ |$$ |  $$ |
# \$$$$$$$\ $$ |\$$$$$$  |$$$$$$$  |  \$$$$  |\$$$$$$$\ $$ |      $$ |$$ |  $$ |\$$$$$$$ |
#  \_______|\__| \______/ \_______/    \____/  \_______|\__|      \__|\__|  \__| \____$$ |
#                                                                               $$\   $$ |
#                                                                               \$$$$$$  |
#                                                                                \______/ 
#____________________________________________________________________________________________


# tag files that contain the extracted barcodes 
files2=list.files(pattern = "*barcodes.txt$")
files2

#___________________________________
# ╔═╗╔═╗╦═╗  ╔═╗╔╗╔╔═╗  ╔═╗╦╦  ╔═╗
# ╠╣ ║ ║╠╦╝  ║ ║║║║║╣   ╠╣ ║║  ║╣ 
# ╚  ╚═╝╩╚═  ╚═╝╝╚╝╚═╝  ╚  ╩╩═╝╚═╝
#__________________________________

library(reclin2, lib.loc="/data/modules/R/4.1.2/lib64/R/library")       # package to minimize memory in R # library(stringdist)    # best package for calculating string distances


dir.create(file.path(here("barcode_clusters", str_replace(files2[1], pattern = "S(\\d+_\\d+).+", replacement = "\\1_clustered"))))

df <- raw.day1      # name the raw.barcodes file into df
df1 <- df %>%           # do not include the barcode length
  dplyr::select(1:2)

df1 <- sample_n(df1, 10000)

my_dist <- function(x, y, ...){          # this is a very cool function where it gives the 1-lv
  1 - stringdist(x, y, method=c("lv"))
}

p <- pair_minsim(df1, on = "barcode", deduplication = TRUE,
                 default_comparator = my_dist, minsim = -2)

group.compmat <- function(p,df1, df){
  g <- graph.data.frame(compmat, directed = TRUE)
  g1 <- compmat %>% 
    mutate(col3 = str_c("group", clusters(g)$membership[as.character(.x)]))
  deduplicate_equivalence(compmat, "group")
  to_add <- data.frame(.x = seq_len(nrow(raw.data)), .y = seq_len(nrow(raw.data)), 
                       simsum = 1)
  
  to_add <- to_add[!(to_add$.x %in% g1$.x), , drop = FALSE]
  to_add <- to_add[!(to_add$.y %in% g1$.y), , drop = FALSE]
  p2 <- bind_rows(g1, to_add)
  attributes(p2) <- attributes(g1)
  
  p3 <- p2 %>% 
    as.data.frame() %>% 
    mutate(
      maxgrp = max(as.integer(gsub("[^0-9]", "", col3)), na.rm = TRUE),
      col3 = if_else(is.na(col3), paste0("group", maxgrp + cumsum(is.na(col3))), col3)
    ) %>%
    select(-maxgrp)
  
  p4 <- p3 %>% 
    mutate(across(1:2, ~raw.data$barcode[match(., raw.data$id)], .names = "barcodes_{.col}")) %>% 
    mutate(across(5:6, .names="{col}.counts", ~ ultra.raw$count[match(., ultra.raw$barcode)])) %>% 
    mutate(lv_distance = case_when(simsum=="0" ~ 1, 
                                   simsum=="-1" ~ 2, 
                                   simsum=="-2" ~ 3,
                                   TRUE ~ 0)) %>% 
    relocate(lv_distance, .after=barcodes_.y) %>% 
    dplyr::rename(groups=col3) %>% 
    dplyr::rename(id=groups,barcode_col1=barcodes_.x,
                  barcode_col2=barcodes_.y,
                  counts_col1=barcodes_.x.counts,
                  counts_col2=barcodes_.y.counts) %>% 
    select(-c(1:3,7)) %>% 
    pivot_longer(cols = -id, names_to = c(".value", "grp"), 
                 names_pattern = "(.*)_(col\\d+)") %>% 
    ungroup() %>% 
    distinct() %>% 
    group_by(id)  %>% 
    mutate(sum_counts = sum(counts)) %>% 
    arrange(id, desc(counts)) %>% 
    mutate(barcodes = paste(barcode, collapse = ",")) %>% 
    dplyr::slice(1) %>% 
    select(id, central_barcode = barcode, barcodes, sum_counts) %>% 
    ungroup() 
  
  
  return(p4)
}

write.table(p4, file=here("barcode_clusters",names(extracts.list[i]),paste0(str_replace(names(extracts.list[i]), pattern = "S(\\d+_\\d+).+", replacement = "\\1_clustered_barcodes"),".txt")))



for(i in 1:lenght(extracts.list)){
  
  dir.create(file.path(here("barcode_clusters", names(extracts.list[i]))))
  
  df <- as.data.frame(extracts.list[i]) %>% 
    dplyr::rename(barcode=1,counts=2,length_barcode=3)
  
  df1 <- df %>% 
    dplyr::select(1:2)
  
  p <- pair_minsim(df1, on = "barcode", deduplication = TRUE,
                   default_comparator = my_dist, minsim = -2)
  
  group.compmat <- function(p,df1, df){
    g <- graph.data.frame(compmat, directed = TRUE)
    g1 <- compmat %>% 
      mutate(col3 = str_c("group", clusters(g)$membership[as.character(.x)]))
    deduplicate_equivalence(compmat, "group")
    to_add <- data.frame(.x = seq_len(nrow(raw.data)), .y = seq_len(nrow(raw.data)), 
                         simsum = 1)
    
    to_add <- to_add[!(to_add$.x %in% g1$.x), , drop = FALSE]
    to_add <- to_add[!(to_add$.y %in% g1$.y), , drop = FALSE]
    p2 <- bind_rows(g1, to_add)
    attributes(p2) <- attributes(g1)
    
    p3 <- p2 %>% 
      as.data.frame() %>% 
      mutate(
        maxgrp = max(as.integer(gsub("[^0-9]", "", col3)), na.rm = TRUE),
        col3 = if_else(is.na(col3), paste0("group", maxgrp + cumsum(is.na(col3))), col3)
      ) %>%
      select(-maxgrp)
    
    p4 <- p3 %>% 
      mutate(across(1:2, ~raw.data$barcode[match(., raw.data$id)], .names = "barcodes_{.col}")) %>% 
      mutate(across(5:6, .names="{col}.counts", ~ ultra.raw$count[match(., ultra.raw$barcode)])) %>% 
      mutate(lv_distance = case_when(simsum=="0" ~ 1, 
                                     simsum=="-1" ~ 2, 
                                     simsum=="-2" ~ 3,
                                     TRUE ~ 0)) %>% 
      relocate(lv_distance, .after=barcodes_.y) %>% 
      dplyr::rename(groups=col3) %>% 
      dplyr::rename(id=groups,barcode_col1=barcodes_.x,
                    barcode_col2=barcodes_.y,
                    counts_col1=barcodes_.x.counts,
                    counts_col2=barcodes_.y.counts) %>% 
      select(-c(1:3,7)) %>% 
      pivot_longer(cols = -id, names_to = c(".value", "grp"), 
                   names_pattern = "(.*)_(col\\d+)") %>% 
      ungroup() %>% 
      distinct() %>% 
      group_by(id)  %>% 
      mutate(sum_counts = sum(counts)) %>% 
      arrange(id, desc(counts)) %>% 
      mutate(barcodes = paste(barcode, collapse = ",")) %>% 
      dplyr::slice(1) %>% 
      select(id, central_barcode = barcode, barcodes, sum_counts) %>% 
      ungroup() 
    
    
    return(p4)
  }
  
  write.table(p4, file=here("barcode_clusters",names(extracts.list[i]),paste0(str_replace(names(extracts.list[i]), pattern = "S(\\d+_\\d+).+", replacement = "\\1_clustered_barcodes"),".txt")))
  
}







