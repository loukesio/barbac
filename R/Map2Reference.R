#' Map and Process Fastq File
#'
#' This function takes a fastq file, maps it to a given reference file,
#' and then processes the result using samtools.
#'
#' @param fastq_file Path to the fastq file.
#' @param reference_file Path to the reference file.
#' @param out_dir Directory to save the processed files.
#' 
#' @return Path to the sorted bam file.
#' @export
#'
#' @examples
#' \dontrun{
#' map_and_process("path/to/fastq", "path/to/reference", "output/directory/")
#' }


map_and_process <- function(fastq_file, reference_file, out_dir) {
  
  # Check if tools are available on system
  if(system("which minimap2", ignore.stdout = TRUE) != 0) {
    stop("Error: 'minimap2' is not found in the system path.")
  }
  
  if(system("which samtools", ignore.stdout = TRUE) != 0) {
    stop("Error: 'samtools' is not found in the system path.")
  }
  
  # Check if files exist
  if(!file.exists(fastq_file)) {
    stop(paste("Error: Fastq file not found at", fastq_file))
  }
  
  if(!file.exists(reference_file)) {
    stop(paste("Error: Reference file not found at", reference_file))
  }
  
  # Check if out directory exists, if not create it
  if(!dir.exists(out_dir)) {
    stop(paste("Error: Output directory not found at", out_dir, ". Please create it or provide a valid directory."))
  }
  
  # Map with minimap2
  sam_file <- paste0(out_dir, basename(fastq_file), ".sam")
  bam_file <- paste0(out_dir, basename(fastq_file), ".bam")
  sorted_bam_file <- paste0(out_dir, basename(fastq_file), "_sorted.bam")
  
  system(sprintf("minimap2 -a %s %s > %s", reference_file, fastq_file, sam_file))
  system(sprintf("samtools view -S -b %s > %s", sam_file, bam_file))
  system(sprintf("samtools sort %s -o %s", bam_file, sorted_bam_file))
  system(sprintf("samtools index %s", sorted_bam_file))
  
  return(sorted_bam_file) # Returning the path to the sorted bam for further usage
}

