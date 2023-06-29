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
library(igraph)
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
    gt::gtsave(filename = here(paste0(str_replace(day1, pattern = "S(\\d+_\\d+).+", replacement = "\\1_summary_rawdata"),".png")))
  
  plot.rawbar <- ggplot(raw.barcodes, aes(x=length_barcode)) +
    geom_histogram(fill="#FF5349", color=NA, alpha=0.8) +
    scale_x_continuous(breaks=seq(0,30, by=5), 
                       limits=c(min(raw.barcodes$length_barcode)-5,max(raw.barcodes$length_barcode) +5)) +
    labs(x="Barcode length (bp)", y="Number of barcodes", 
         title=paste0(str_replace(day1, pattern = "S(\\d+_\\d+).+", replacement = "\\1_raw_data")), caption=c("Barbac: A versatile tool to quantify barcodes")) + # find a way to make it more general
    theme_bw() 
  
  ggsave(here(paste0(str_replace(day1, pattern = "S(\\d+_\\d+).+", replacement = "\\1_barcode_length"),".png")))
  
  
  return(raw.barcodes)
}
raw.day1 <- barcodes.extract.m1(day1)
raw.day1

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

dir.create(file.path("/home/theodosiou/Projects/Barcodes/Static_Experiment/Amplicon_Barcodes/rep2/merged/assembled_reads/minimap2_output/barcode_clusters", str_replace(files2[1], pattern = "S(\\d+_\\d+).+", replacement = "\\1_clustered")))
setwd(file.path("/home/theodosiou/Projects/Barcodes/Static_Experiment/Amplicon_Barcodes/rep2/merged/assembled_reads/minimap2_output/barcode_clusters",str_replace(files2[1], pattern = "S(\\d+_\\d+).+", replacement = "\\1_clustered")))

df <- raw.day1 %>% 
  mutate(id=row_number())                # name the raw.barcodes file into df

df1 <- df %>%                            # I am taking a subsample with only the barcodes column and the first 11000 rows.
  dplyr::select(1) 

my_dist <- function(x, y, ...){          # this is a very cool function where it gives the 1-lv
  1 - stringdist(x, y, method=c("lv"))
}

p <- pair_minsim(df1, on = "barcode", deduplication = TRUE,
                 default_comparator = my_dist, minsim = -2)

#__________________________________________________________________
#   __ _  _ __  ___   _   _  _ __     ___  ___   _ __ ___   _ __  
#  / _` || '__|/ _ \ | | | || '_ \   / __|/ _ \ | '_ ` _ \ | '_ \ 
# | (_| || |  | (_) || |_| || |_) |_| (__| (_) || | | | | || |_) |
#  \__, ||_|   \___/  \__,_|| .__/(_)\___|\___/ |_| |_| |_|| .__/ 
#   __/ |                   | |                            | |    
#  |___/                    |_|                            |_|    
#__________________________________________________________________

group.compmat <- function(p,df1, df){
  
  g <- graph.data.frame(p, directed = TRUE)         # this function finds the indirect connections
  g1 <- p %>% 
    mutate(col3 = str_c("group", clusters(g)$membership[as.character(.x)])) # g1 allows to connect indirect connections to the same group
  
  deduplicate_equivalence(p, "group")                                       # this function adds the barcodes that did not come out of the comparison matrix
  to_add <- data.frame(.x = seq_len(nrow(df1)), .y = seq_len(nrow(df1)),    # here I am adding the id and the dinstance of these barcodes which is zero
                       simsum = 1)
  to_add <- to_add[!(to_add$.x %in% g1$.x), , drop = FALSE]                 # need to find out what this one does, but its so cool !!!!!
  to_add <- to_add[!(to_add$.y %in% g1$.y), , drop = FALSE]
  
  p2 <- bind_rows(g1, to_add)
  attributes(p2) <- attributes(g1)
  
  
  p3 <- p2 %>%                                                              # you make p2 to data.frame and you add the group word inside
    as.data.frame() %>% 
    mutate(
      maxgrp = max(as.integer(gsub("[^0-9]", "", col3)), na.rm = TRUE),
      col3 = if_else(is.na(col3), paste0("group", maxgrp + cumsum(is.na(col3))), col3)
    ) %>%
    select(-maxgrp)
  
  
  p4 <- p3 %>%                        
    mutate(across(1:2, ~df$barcode[match(., df$id)], .names = "barcodes_{.col}")) %>% 
    mutate(across(5:6, .names="{col}.counts", ~ df$counts[match(., df$barcode)])) %>% 
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
    group_by(id) %>% 
    distinct(barcode, .keep_all = TRUE) %>% 
    mutate(sum_counts = sum(counts)) %>% 
    mutate(barcodes = paste(barcode, collapse = ",")) %>% 
    select(id, central_barcode = barcode, barcodes, sum_counts, counts) %>% 
    ungroup() 
  
  return(p4)
}

barcode <- group.compmat(p,df1,df)

write.table(barcode, file=paste0(str_replace(files2[1], pattern = "S(\\d+_\\d+).+", replacement = "\\1_barcodes"),".txt"))


