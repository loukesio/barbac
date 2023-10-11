#' Extract Barcodes from BAM File
#'
#' This function extracts barcode sequences from a BAM file based on the specified reference name and position range.
#'
#' @param bam_file Character string. Path to the BAM file from which to extract barcodes.
#' @param ref_name Character string. Name of the reference to be used. Default is "Reference_barcodes".
#' @param start_pos Numeric. Start position of the barcode in the reference. Default is 54.
#' @param end_pos Numeric. End position of the barcode in the reference. Default is 78.
#' @return Character string. Path to the generated barcode file.
#' @import dplyr
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsamtools stackStringsFromBam
#' @importFrom stringr str_length
#' @export
extract_barcodes <- function(bam_file, ref_name="Reference_barcodes", start_pos=54, end_pos=78) {
  
  # Check if bam file exists
  if (!file.exists(bam_file)) {
    stop(paste("Error: BAM file not found at", bam_file))
  }
  
  # Ensure that start_pos and end_pos are valid
  if (!is.numeric(start_pos) || !is.numeric(end_pos) || start_pos >= end_pos) {
    stop("Error: Invalid start or end position values.")
  }
  
  param_range <- sprintf("%s:%d-%d", ref_name, start_pos, end_pos)
  
  stack <- stackStringsFromBam(bam_file,
                               param=GRanges(param_range),
                               use.names = TRUE,
                               what = "seq")
  
  data <- stack %>%
    as.data.frame() %>%
    dplyr::rename(barcode=x) %>%
    group_by(barcode) %>%
    dplyr::count() %>%
    arrange(desc(n)) %>%
    dplyr::rename(counts=n) %>%
    mutate(length_barcode=str_length(barcode))
  
  # Generate the output file name
  output_file <- paste0(tools::file_path_sans_ext(bam_file), "_barcodes.txt")
  write.table(data, output_file, quote=FALSE, sep=",", row.names=FALSE)
  
  return(output_file) # Return the path to the extracted barcode file for user's reference
}
