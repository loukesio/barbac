#' Extract Barcodes from BAM File
#'
#' This function extracts barcode sequences from a BAM file based on the specified 
#' reference name and position range, and saves the result as a clean CSV file.
#'
#' @param bam_file Character string. Path to the BAM file from which to extract barcodes.
#' @param ref_name Character string. Name of the reference to be used. Default is "Reference_barcodes".
#' @param start_pos Numeric. Start position of the barcode in the reference. Default is 54.
#' @param end_pos Numeric. End position of the barcode in the reference. Default is 78.
#' @param output_file Character string or NULL. Path for output CSV. If NULL, creates filename 
#'   based on input BAM file. Default is NULL.
#' @param min_count Numeric. Minimum count threshold to include a barcode. Default is 1.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#'
#' @return A character string with the path to the generated CSV file.
#' 
#' @details
#' This function:
#' \itemize{
#'   \item Reads aligned sequences from a BAM file at specified genomic coordinates
#'   \item Extracts barcode sequences from the specified region
#'   \item Counts occurrence of each unique barcode
#'   \item Filters by minimum count threshold
#'   \item Saves results to CSV with columns: barcode, counts, barcode_length
#' }
#'
#' @examples
#' \dontrun{
#' # Extract barcodes from default region
#' barbac_xtr("sample1_sorted.bam")
#' 
#' # Extract from custom region
#' barbac_xtr(
#'   bam_file = "sample1_sorted.bam",
#'   ref_name = "chr1",
#'   start_pos = 100,
#'   end_pos = 125
#' )
#' 
#' # Specify output file and filter low counts
#' barbac_xtr(
#'   bam_file = "sample1_sorted.bam",
#'   output_file = "sample1_filtered_barcodes.csv",
#'   min_count = 5
#' )
#' }
#'
#' @import dplyr
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicAlignments stackStringsFromBam
#' @importFrom stringr str_length
#' @importFrom readr write_csv
#' @export
barbac_xtr <- function(bam_file, 
                       ref_name = "Reference_barcodes", 
                       start_pos = 54, 
                       end_pos = 78,
                       output_file = NULL,
                       min_count = 1,
                       verbose = TRUE) {
  
  # -----------------------------
  # Input validation
  # -----------------------------
  if (!file.exists(bam_file)) {
    stop("\u274C BAM file not found at: ", bam_file)
  }
  
  # Check for BAM index
  bai_file <- paste0(bam_file, ".bai")
  if (!file.exists(bai_file)) {
    warning("\u26A0 BAM index (.bai) not found. Creating index...")
    system2("samtools", c("index", shQuote(bam_file)))
  }
  
  if (!is.numeric(start_pos) || !is.numeric(end_pos) || start_pos >= end_pos) {
    stop("\u274C Invalid start or end position values. start_pos must be < end_pos.")
  }
  
  if (!is.numeric(min_count) || min_count < 1) {
    stop("\u274C min_count must be a positive integer.")
  }
  
  # -----------------------------
  # Define output file
  # -----------------------------
  if (is.null(output_file)) {
    output_file <- paste0(tools::file_path_sans_ext(bam_file), "_barcodes.csv")
  }
  
  # Ensure output directory exists
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # -----------------------------
  # Extract barcodes
  # -----------------------------
  if (verbose) {
    message("\u25B6 Extracting barcodes from BAM file...")
    message(sprintf("  BAM file: %s", basename(bam_file)))
    message(sprintf("  Reference: %s", ref_name))
    message(sprintf("  Region: %d-%d (length: %d bp)", start_pos, end_pos, end_pos - start_pos + 1))
  }
  
  # Define genomic range for extraction
  param_range <- sprintf("%s:%d-%d", ref_name, start_pos, end_pos)
  
  # Extract sequences from BAM
  tryCatch({
    stack <- GenomicAlignments::stackStringsFromBam(
      bam_file,
      param = GenomicRanges::GRanges(param_range),
      use.names = TRUE,
      what = "seq"
    )
  }, error = function(e) {
    stop("\u274C Failed to extract sequences from BAM file.\n",
         "  Error: ", e$message, "\n",
         "  Check that reference name '", ref_name, "' exists in BAM file.")
  })
  
  if (length(stack) == 0) {
    warning("\u26A0 No sequences found in the specified region!")
    # Create empty output
    data <- tibble::tibble(
      barcode = character(),
      counts = integer(),
      barcode_length = integer()
    )
  } else {
    # -----------------------------
    # Process and summarize barcodes
    # -----------------------------
    if (verbose) message("\u25B6 Processing barcode sequences...")
    
    data <- stack %>%
      as.data.frame() %>%
      dplyr::rename(barcode = x) %>%
      dplyr::mutate(barcode = as.character(barcode)) %>%
      dplyr::group_by(barcode) %>%
      dplyr::summarise(counts = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(counts)) %>%
      dplyr::mutate(barcode_length = stringr::str_length(barcode))
    
    # Apply count filter
    if (min_count > 1) {
      n_before <- nrow(data)
      data <- data %>% dplyr::filter(counts >= min_count)
      n_after <- nrow(data)
      
      if (verbose && n_before > n_after) {
        message(sprintf("  Filtered out %d barcodes with counts < %d", 
                        n_before - n_after, min_count))
      }
    }
  }
  
  # -----------------------------
  # Save results
  # -----------------------------
  readr::write_csv(data, output_file)
  
  if (verbose) {
    message("\u2705 Barcode extraction complete!")
    message(sprintf("  Total unique barcodes: %s", format(nrow(data), big.mark = ",")))
    if (nrow(data) > 0) {
      message(sprintf("  Total reads: %s", format(sum(data$counts), big.mark = ",")))
      message(sprintf("  Most abundant barcode: %s (count: %s)", 
                      data$barcode[1], format(data$counts[1], big.mark = ",")))
      message(sprintf("  Barcode length range: %d - %d bp", 
                      min(data$barcode_length), max(data$barcode_length)))
    }
    message(sprintf("  Output saved to: %s", output_file))
  }
  
  return(invisible(output_file))
}


#' Generate Summary Statistics and Plots for Extracted Barcodes
#'
#' @param file A data frame or path to CSV file with columns: barcode, counts, barcode_length
#' @param barcode_length A numeric vector of length 2 specifying desired barcode length range [min, max]
#' @param fill_color Histogram fill color. Default is "#FF5349" (red).
#' @param save_plot Logical or character. If TRUE, saves to default name. If character, saves to that path. Default is FALSE.
#' @param plot_width Numeric. Width of saved plot in inches. Default is 12.
#' @param plot_height Numeric. Height of saved plot in inches. Default is 8.
#' @param verbose Logical. Print summary statistics. Default is TRUE.
#' 
#' @return Patchwork plot object with histograms and summary table
#'
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @importFrom gridExtra tableGrob
#' @importFrom readr read_csv
#' @importFrom stats median
#' @export
barbac_xtr.stats <- function(file,
                             barcode_length, 
                             fill_color = "#FF5349",
                             save_plot = FALSE,
                             plot_width = 12,
                             plot_height = 8,
                             verbose = TRUE) {
  
  # Load patchwork if not already loaded
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install it with: install.packages('patchwork')")
  }
  
  # -----------------------------
  # Input validation
  # -----------------------------
  if (!is.numeric(barcode_length) || length(barcode_length) != 2) {
    stop("\u274C 'barcode_length' should be a numeric vector of length 2 (e.g., c(20, 30)).")
  }
  
  if (barcode_length[1] >= barcode_length[2]) {
    stop("\u274C First value in 'barcode_length' must be less than the second.")
  }
  
  # -----------------------------
  # Read data
  # -----------------------------
  if (is.character(file)) {
    if (!file.exists(file)) {
      stop("\u274C File not found: ", file)
    }
    if (verbose) message("\u25B6 Reading barcode data from: ", basename(file))
    data <- readr::read_csv(file, show_col_types = FALSE)
    input_file <- file
  } else {
    data <- file
    input_file <- NULL
  }
  
  # Validate required columns
  required <- c("barcode", "counts", "barcode_length")
  if (!all(required %in% colnames(data))) {
    stop("\u274C Input must contain columns: ", paste(required, collapse = ", "))
  }
  
  # -----------------------------
  # Print summary statistics
  # -----------------------------
  if (verbose) {
    message("\n=== BARCODE SUMMARY STATISTICS ===")
    message(sprintf("Total unique barcodes: %s", format(nrow(data), big.mark = ",")))
    message(sprintf("Total reads: %s", format(sum(data$counts), big.mark = ",")))
    message(sprintf("Barcode length range: %d - %d bp", 
                    min(data$barcode_length), max(data$barcode_length)))
    message(sprintf("Mean barcode length: %.1f bp", mean(data$barcode_length)))
    message(sprintf("Median count: %s", format(median(data$counts), big.mark = ",")))
    message(sprintf("Mean count: %.1f", mean(data$counts)))
  }
  
  # -----------------------------
  # Create bin summary table
  # -----------------------------
  bin_data <- data %>%
    dplyr::mutate(
      bin = dplyr::case_when(
        barcode_length < barcode_length[1] ~ 
          sprintf("< %d bp", barcode_length[1]),
        barcode_length >= barcode_length[1] & barcode_length <= barcode_length[2] ~ 
          sprintf("%d - %d bp", barcode_length[1], barcode_length[2]),
        TRUE ~ 
          sprintf("> %d bp", barcode_length[2])
      )
    ) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      barcodes = dplyr::n(),
      reads = sum(counts),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      `% barcodes` = sprintf("%.1f%%", 100 * barcodes / sum(barcodes)),
      `% reads` = sprintf("%.1f%%", 100 * reads / sum(reads))
    )
  
  # Format for table
  bin_table <- bin_data %>%
    dplyr::mutate(
      barcodes = format(barcodes, big.mark = ","),
      reads = format(reads, big.mark = ",")
    ) %>%
    dplyr::rename(
      `Length bin` = bin,
      `Barcodes` = barcodes,
      `Reads` = reads
    )
  
  t1 <- gridExtra::tableGrob(bin_table, rows = NULL)
  
  if (verbose) {
    message("\n=== LENGTH DISTRIBUTION ===")
    print(bin_data, n = Inf)
  }
  
  # -----------------------------
  # Plot 1: Barcode length distribution
  # -----------------------------
  p1 <- ggplot2::ggplot(data, ggplot2::aes(x = barcode_length)) +
    ggplot2::geom_histogram(fill = fill_color, color = NA, alpha = 0.8, bins = 30) +
    ggplot2::geom_vline(xintercept = barcode_length[1], 
                        linetype = "dashed", color = "gray30", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = barcode_length[2], 
                        linetype = "dashed", color = "gray30", linewidth = 0.5) +
    ggplot2::labs(
      x = "Barcode length (bp)", 
      y = "Number of barcodes", 
      title = "Barcode Length Distribution"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # -----------------------------
  # Plot 2: Barcode count distribution (log scale)
  # -----------------------------
  p2 <- ggplot2::ggplot(data, ggplot2::aes(x = log10(counts))) +
    ggplot2::geom_histogram(fill = fill_color, color = NA, alpha = 0.8, bins = 30) +
    ggplot2::labs(
      x = "log10(Barcode counts)", 
      y = "Number of barcodes", 
      title = "Barcode Count Distribution"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # -----------------------------
  # Plot 3: Shannon entropy distribution
  # -----------------------------
  shannon_entropy <- function(seq) {
    n <- nchar(seq)
    if (n == 0) return(0)
    nucleotides <- strsplit(seq, "")[[1]]
    probs <- table(nucleotides) / n
    -sum(probs * log2(probs))
  }
  
  if (verbose) message("\n\u25B6 Calculating Shannon entropy...")
  
  data_with_entropy <- data %>%
    dplyr::mutate(entropy = sapply(barcode, shannon_entropy))
  
  p3 <- ggplot2::ggplot(data_with_entropy, ggplot2::aes(x = entropy)) +
    ggplot2::geom_histogram(fill = fill_color, color = NA, alpha = 0.8, bins = 30) +
    ggplot2::labs(
      x = "Shannon Entropy", 
      y = "Number of barcodes", 
      title = "Barcode Entropy Distribution"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (verbose) {
    message(sprintf("Mean entropy: %.3f", mean(data_with_entropy$entropy)))
    message(sprintf("Entropy range: %.3f - %.3f", 
                    min(data_with_entropy$entropy), 
                    max(data_with_entropy$entropy)))
  }
  
  # -----------------------------
  # Combine plots using patchwork
  # -----------------------------
  combined_plot <- (p1 | p2) / (p3 | patchwork::wrap_elements(t1))
  
  # -----------------------------
  # Save plot if requested
  # -----------------------------
  if (!isFALSE(save_plot)) {
    if (isTRUE(save_plot)) {
      # Generate default filename
      if (!is.null(input_file)) {
        plot_file <- paste0(tools::file_path_sans_ext(input_file), "_stats.png")
      } else {
        plot_file <- "barcode_stats.png"
      }
    } else {
      plot_file <- save_plot
    }
    
    ggplot2::ggsave(
      filename = plot_file,
      plot = combined_plot,
      width = plot_width,
      height = plot_height,
      dpi = 300
    )
    
    if (verbose) message(sprintf("\n\u2705 Plot saved to: %s", plot_file))
  }
  
  if (verbose) message("\n\u2705 Analysis complete!\n")
  
  return(combined_plot)
}

