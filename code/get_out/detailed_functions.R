table <- function(file, barcode_length) {
  # Check if the 'barcode_length' parameter is a numeric vector of length 2
  if (!is.numeric(barcode_length) || length(barcode_length) != 2) {
    stop("Error: 'barcode_length' should be a numeric vector of length 2.")
  }
  
  # Check if 'file' is a data frame or tibble
  if (!is.data.frame(file) && !is.tbl(file)) {
    stop("Error: 'file' should be a data frame or tibble.")
  }

  file %>%
    mutate(
      bin = case_when(
        (length_barcode < barcode_length[1] & length_barcode > 0) ~ paste0("(", 0, ",", barcode_length[1], ")"),
        (length_barcode >= barcode_length[1] & length_barcode <= barcode_length[2]) ~ paste0("[", barcode_length[1], ",", barcode_length[2], "]"),
        TRUE ~ paste0("[", barcode_length[2] + 1, ",∞)")
      )
    ) %>%
    group_by(bin) %>%
    count() %>%
    gridExtra::tableGrob()
}


t1 <- table(file, barcode_length)

############################
# plot 1 -- barcode length
############################

p1 <- 
  test %>% 
  ggplot(., aes(x=length_barcode)) +
  geom_histogram(fill="#FF5349", color=NA, alpha=0.8) +
  labs(x="Barcode length (bp)", y="Number of barcodes", title="Barcode length distribution") + 
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

############################
# plot 2 -- barcode counts
############################

p2 <- 
  test %>%   
  ggplot(., aes(x=log(counts))) +
  geom_histogram(fill="#FF5349", color=NA, alpha=0.8) +
  labs(x="log(Barcode counts)", y="Number of barcodes", title="Barcode counts distribution") + 
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

##############################
# plot 3 -- sequence entropy
##############################

# Function to calculate Shannon entropy
shannon_entropy <- function(seq) {
  n <- nchar(seq)
  nucleotides <- unique(strsplit(seq, "")[[1]])
  counts <- sapply(nucleotides, function(nuc) sum(strsplit(seq, "")[[1]] == nuc))
  probabilities <- counts / n
  entropy <- -sum(probabilities * log2(probabilities))
  return(entropy)
}

p3 <- 
  test %>% 
  mutate(entropy = sapply(barcode, shannon_entropy)) %>% 
  ggplot() +
  geom_histogram(aes(x=entropy), fill="#FF5349") +
  theme_bw()


(p1 | p2) / (p3 | t1)
