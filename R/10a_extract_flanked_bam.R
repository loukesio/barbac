# Extract actual query bases, with reference coverage providing locus selection.
# No CIGAR insertion/deletion is rewritten into a fixed-width barcode string.
.barbac_xtr_flanks <- function(bam_file, ref_name, start_pos, end_pos,
                               output_file, min_count, verbose, flank_pattern,
                               barcode_group, read_window, reverse_complement,
                               include_read_ids, yield_size) {
  scalar_integer <- function(x, minimum) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x >= minimum && x == floor(x)
  if (!file.exists(bam_file)) stop("BAM file not found: ", bam_file)
  if (!scalar_integer(start_pos, 2) || !scalar_integer(end_pos, start_pos) ||
      !scalar_integer(min_count, 1) || !scalar_integer(barcode_group, 1) ||
      !scalar_integer(yield_size, 1)) stop("Invalid flank-extraction coordinates or integer settings.")
  if (length(ref_name) != 1L || is.na(ref_name) || !nzchar(ref_name) ||
      length(flank_pattern) != 1L || is.na(flank_pattern) || !nzchar(flank_pattern)) {
    stop("Supply one reference name and one nonempty flank pattern.")
  }
  for (x in list(reverse_complement, include_read_ids, verbose)) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) stop("Logical settings must be TRUE or FALSE.")
  }
  if (!is.null(read_window) && (length(read_window) != 2L ||
      !scalar_integer(read_window[1], 1) || !scalar_integer(read_window[2], read_window[1]))) {
    stop("read_window must contain two positive, increasing inclusive query coordinates.")
  }
  if (include_read_ids && min_count != 1) stop("Read-level extraction requires min_count = 1.")
  probe <- regexpr(flank_pattern, "", perl = TRUE)
  if (is.null(attr(probe, "capture.start")) || ncol(attr(probe, "capture.start")) < barcode_group) {
    stop("barcode_group does not exist in flank_pattern.")
  }
  targets <- Rsamtools::scanBamHeader(bam_file)[[1]]$targets
  if (!ref_name %in% names(targets)) stop("Reference not found in BAM: ", ref_name)
  if (end_pos >= targets[[ref_name]]) stop("Barcode region requires a reference base on either side.")
  if (is.null(output_file)) output_file <- paste0(tools::file_path_sans_ext(bam_file), "_barcodes.csv")
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("barbac-extract-", tmpdir = dirname(output_file))
  on.exit(unlink(temporary), add = TRUE)
  empty_reads <- data.frame(read_id = character(), barcode = character(), barcode_length = integer())
  if (include_read_ids) readr::write_csv(empty_reads, temporary)
  bam <- Rsamtools::BamFile(bam_file, yieldSize = as.integer(yield_size))
  open(bam)
  on.exit(close(bam), add = TRUE)
  param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE),
    what = c("qname", "rname", "pos", "cigar", "seq"), reverseComplement = FALSE)
  stats <- c(primary_mapped_alignments = 0, covering_locus = 0,
             pattern_failed = 0, matched_alignments = 0)
  chunks <- list()
  repeat {
    part <- Rsamtools::scanBam(bam, param = param)[[1]]
    if (!length(part$qname)) break
    stats[[1]] <- stats[[1]] + length(part$qname)
    cigars <- unique(part$cigar)
    widths <- vapply(regmatches(cigars, gregexpr("[0-9]+[MDN=X]", cigars)), function(ops) {
      sum(as.numeric(sub("[MDN=X]$", "", ops)))
    }, numeric(1))
    right <- part$pos + widths[match(part$cigar, cigars)] - 1L
    keep <- as.character(part$rname) == ref_name & part$pos <= start_pos - 1L & right >= end_pos + 1L &
      !grepl("[HN]", part$cigar)
    stats[[2]] <- stats[[2]] + sum(keep)
    if (!any(keep)) next
    sequences <- part$seq[keep]
    if (reverse_complement) sequences <- Biostrings::reverseComplement(sequences)
    sequences <- as.character(sequences)
    if (!is.null(read_window)) sequences <- substr(sequences, read_window[1], read_window[2])
    hits <- regexpr(flank_pattern, sequences, perl = TRUE)
    positions <- attr(hits, "capture.start")[, barcode_group]
    lengths <- attr(hits, "capture.length")[, barcode_group]
    matched <- positions > 0L & lengths > 0L
    stats[[3]] <- stats[[3]] + sum(!matched)
    stats[[4]] <- stats[[4]] + sum(matched)
    if (!any(matched)) next
    barcodes <- substring(sequences[matched], positions[matched], positions[matched] + lengths[matched] - 1L)
    if (include_read_ids) {
      readr::write_csv(data.frame(read_id = part$qname[keep][matched], barcode = barcodes,
                                 barcode_length = nchar(barcodes)), temporary, append = TRUE)
    } else chunks[[length(chunks) + 1L]] <- table(barcodes)
  }
  if (!include_read_ids) {
    if (length(chunks)) {
      values <- unlist(chunks)
      totals <- tapply(as.numeric(values), names(values), sum)
      data <- data.frame(barcode = names(totals), counts = as.numeric(totals))
      data <- data[data$counts >= min_count, , drop = FALSE]
      data <- data[order(-data$counts, data$barcode), , drop = FALSE]
      data$barcode_length <- nchar(data$barcode)
    } else data <- data.frame(barcode = character(), counts = numeric(), barcode_length = integer())
    readr::write_csv(data, temporary)
  }
  if (!file.rename(temporary, output_file)) stop("Could not finalize extraction output: ", output_file)
  if (verbose) message("Extracted ", stats[[4]], " barcode observations from ", stats[[2]], " covering alignments.")
  structure(output_file, extraction_stats = as.list(stats))
}
