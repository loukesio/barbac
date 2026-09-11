# Read both observed barcode components from primary merged BAM alignments.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
if (nzchar(Sys.getenv("BARBAC_R_LIBRARY"))) .libPaths(c(Sys.getenv("BARBAC_R_LIBRARY"), .libPaths()))
library(barbac)
if (!"flank_pattern" %in% names(formals(barbac_xtr))) stop("Install the current branch of barbac before BAM extraction")
bam <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out, "bam_components.json"))) stop("Completed component extraction exists")
pattern <- paste0("\\D*?(GTACC|GGACC|GGTCC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGTAC.)",
  "(\\D{24,28})(.TAACT|A.AACT|AT.ACT|ATA.CT|ATAA.T|ATAAC|AAACT|ATACT|ATAAT)\\D*")
started <- proc.time()[["elapsed"]]
d <- barbac_xtr(bam, "chen2023_barcode_amplicon", 50, 75, file.path(out, "diverse_reads.csv"),
  flank_pattern = pattern, barcode_group = 2L, read_window = c(54, 99), include_read_ids = TRUE)
e <- barbac_xtr(bam, "chen2023_barcode_amplicon", 110, 135, file.path(out, "environment_reads.csv"),
  flank_pattern = pattern, barcode_group = 2L, read_window = c(40, 85),
  reverse_complement = TRUE, include_read_ids = TRUE)
jsonlite::write_json(list(status = "complete", seconds = proc.time()[["elapsed"]] - started,
  diverse = attr(d, "extraction_stats"), environment = attr(e, "extraction_stats"),
  source = "mapped PEAR consensus query sequences; flanks selected inside the published query windows",
  barcode_group = 2L, pattern = pattern, library_paths = .libPaths(), bam_md5 = unname(tools::md5sum(bam)),
  output_md5 = as.list(tools::md5sum(c(d, e)))),
  file.path(out, "bam_components.json"), pretty = TRUE, auto_unbox = TRUE)
