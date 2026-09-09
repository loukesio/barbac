# Run the existing barbac pipeline on the first complete, verified sample.
# Usage: Rscript run_mapping_pilot.R /absolute/path/to/chen_work
# This stops at BAM generation; it does not claim paired/UMI barcode extraction.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply one work directory containing raw/SRR22757105/")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
kit <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
work <- normalizePath(args[[1]], mustWork = TRUE)
manifest <- read.delim(file.path(kit, "samples.tsv"), check.names = FALSE)
sample <- manifest[1L, , drop = FALSE]
raw <- file.path(work, "raw", sample$run)
files <- file.path(raw, paste0(sample$run, "_", 1:2, ".fastq.gz"))
if (!all(file.exists(files))) stop("Download the first sample with run_sample.py download first")
expected_md5 <- unlist(sample[c("r1_md5", "r2_md5")], use.names = FALSE)
expected_bytes <- as.numeric(unlist(sample[c("r1_bytes", "r2_bytes")], use.names = FALSE))
if (!identical(unname(tools::md5sum(files)), expected_md5) ||
    !all(file.info(files)$size == expected_bytes)) stop("FASTQ size/MD5 mismatch")
reference <- file.path(kit, "reference", "chen2023_masked_amplicon.fasta")
if (!file.exists(reference)) stop("Masked candidate reference is missing")
output <- file.path(work, "full_pipeline")
if (dir.exists(output)) stop("Output already exists; use a new work directory for a rerun")

library(barbac)
use_barbac_env()
check_barbac_tools()
started <- proc.time()[["elapsed"]]
result <- run_cli_pipeline(
  data.frame(sample = sample$run, R1 = files[[1]], R2 = files[[2]]),
  reference = reference, output_dir = output, verbose = FALSE)
seconds <- proc.time()[["elapsed"]] - started
bam <- file.path(output, "merged", "bam",
                 paste0(sample$run, "_ANC.assembled_sorted.bam"))
expected <- c(bam, paste0(bam, ".bai"), file.path(output, "bam_summary.csv"),
              file.path(output, "fastQC", paste0(sample$run, "_", 1:2, "_fastqc.zip")))
if (!all(file.exists(expected)) || any(file.info(expected)$size == 0)) {
  stop("Pipeline did not produce every expected output; inspect pipeline.log")
}
if (system2(Sys.which("samtools"), c("quickcheck", shQuote(bam))) != 0L) {
  stop("BAM integrity check failed")
}
if (any(result$stats$mapped <= 0)) stop("No mapped reads; inspect the pipeline outputs")
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines(capture.output(dput(list(
  status = "BAM generation complete; paired/UMI extraction not run",
  pipeline_seconds = seconds, sample = sample, stats = result$stats,
  reference = reference, reference_md5 = tools::md5sum(reference),
  tools = Sys.which(c("fastqc", "pear", "minimap2", "samtools", "multiqc")),
  barbac_version = as.character(packageVersion("barbac"))))),
  file.path(output, "mapping_receipt.R"))
print(result$stats)
message(sprintf("BAM pipeline completed in %.2f seconds; extraction and clustering remain.", seconds))
