#' Summarize Barcode Statistics
#'
#' This function generates a summary of barcode statistics including barcode length distribution,
#' barcode counts distribution, and barcode entropy distribution.
#'
#' @param file A data frame or tibble containing the barcode data.
#' @param barcode_length A numeric vector of length 2 specifying the minimum and maximum barcode lengths.
#' @param fill_color The fill color for the histograms.
#'
#' @return A combined plot of the barcode length distribution, barcode counts distribution,
#' and barcode entropy distribution.
#'
#' @import dplyr
#' @import ggplot2
#' @import gridExtra
#'
#' @keywords barcode statistics plot histogram entropy
#'
#' @examples
#' # Generate summary statistics for barcode data
#' barbac_xtr.stats(my_data, c(5, 15))
#'
#' @export
#'
barbac_xtr.stats <- function(file, barcode_length, fill_color = "#FF5349") {
  # Check if the barcode_length parameter is a numeric vector of length 2
  if (!is.numeric(barcode_length) || length(barcode_length) != 2) {
    stop("Error: 'barcode_length' should be a numeric vector of length 2.")
  }
  
  # Check if 'file' is a data frame or tibble
  if (!is.data.frame(file) && !is.data.frame(file)) {
    stop("Error: 'file' should be a data frame or tibble.")
  }
  
  t1 <- 
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
    dplyr::rename(barcodes=n) %>% 
    gridExtra::tableGrob(., rows = NULL, theme = ttheme_default(core = list(bg_params = list(fill = "white")), colhead = list(fg_params = list(col = "black"), bg_params = list(fill = "white", lwd = 1, lty = 1, lineend = "butt", col = "#333333"))))
  
  p1 <- ggplot2::ggplot(file, aes(x=length_barcode)) +
    ggplot2::geom_histogram(fill = fill_color, color = NA, alpha = 0.8) +
    ggplot2::labs(x="Barcode length (bp)", y="Number of barcodes", title="Barcode length distribution") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5))
  
  p2 <- ggplot2::ggplot(file, aes(x=log(counts))) +
    ggplot2::geom_histogram(fill = fill_color, color = NA, alpha = 0.8) +
    ggplot2::labs(x="log(Barcode counts)", y="Number of barcodes", title="Barcode counts distribution") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5))
  
  shannon_entropy <- function(seq) {
    n <- nchar(seq)
    nucleotides <- unique(strsplit(seq, "")[[1]])
    counts <- sapply(nucleotides, function(nuc) sum(strsplit(seq, "")[[1]] == nuc))
    probabilities <- counts / n
    entropy <- -sum(probabilities * log2(probabilities))
    return(entropy)
  }
  
  p3 <- dplyr::mutate(file, entropy = sapply(barcode, shannon_entropy)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(aes(x=entropy), fill = fill_color) +
    ggplot2::theme_bw()
  
  return((p1 | p2) / (p3 | t1))
}

## Quiet concerns about global variables
utils::globalVariables(c("%>%", "aes", "barcode", "bin", "counts", "is.tbl", "length_barcode", "n","entropy"))
