# Existing barbac FASTQ -> BAM workflow for a checked manifest row.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
kit <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
work <- normalizePath(args[1], mustWork = TRUE)
index <- as.integer(args[2]) + 1L
rows <- read.delim(file.path(kit, "samples.tsv"), check.names = FALSE)
stopifnot(index >= 1L, index <= nrow(rows))
row <- rows[index, , drop = FALSE]
files <- file.path(work, "raw", row$run, paste0(row$run, "_", 1:2, ".fastq.gz"))
stopifnot(all(file.exists(files)))
stopifnot(identical(unname(tools::md5sum(files)), unlist(row[c("r1_md5", "r2_md5")], use.names = FALSE)))
stopifnot(all(file.info(files)$size == as.numeric(unlist(row[c("r1_bytes", "r2_bytes")]))))
if (nzchar(Sys.getenv("BARBAC_R_LIBRARY"))) .libPaths(c(Sys.getenv("BARBAC_R_LIBRARY"), .libPaths()))
library(barbac)
use_barbac_env()
output <- file.path(work, "mapping", row$sample)
if (dir.exists(output)) stop("Mapping output exists; verify its receipt or use a new directory")
reference <- file.path(kit, "reference", "chen2023_masked_amplicon.fasta")
started <- proc.time()[["elapsed"]]
result <- run_cli_pipeline(data.frame(sample = row$run, R1 = files[1], R2 = files[2]),
  reference = reference, output_dir = output, verbose = FALSE)
bam <- file.path(output, "merged", "bam", paste0(row$run, "_ANC.assembled_sorted.bam"))
required <- c(bam, paste0(bam, ".bai"), file.path(output, "fastQC", paste0(row$run, "_", 1:2, "_fastqc.zip")))
stopifnot(all(file.exists(required)), all(file.info(required)$size > 0),
  system2(Sys.which("samtools"), c("quickcheck", shQuote(bam))) == 0L,
  all(result$stats$mapped > 0))
jsonlite::write_json(list(status = "complete", stage = "mapping", sample = row,
  bam = bam, reference_md5 = unname(tools::md5sum(reference)), input_md5 = unname(tools::md5sum(files)),
  seconds = proc.time()[["elapsed"]] - started, stats = result$stats,
  commands = result$commands, barbac_version = as.character(packageVersion("barbac"))),
  file.path(output, "mapping.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
print(result$stats)
