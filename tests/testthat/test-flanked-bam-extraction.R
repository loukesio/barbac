flanked_bam_fixture <- function() {
  samtools <- .barbac_tool("samtools")
  skip_if(!nzchar(samtools), "samtools unavailable")
  folder <- tempfile("barbac-flank-test-")
  dir.create(folder)
  prefix <- "ACGTCAGGTACC"
  suffix <- "ATAACTGGTACG"
  barcode <- "ACGTACGTACGTACGTACGTACGTAC"
  bcs <- c(barcode, paste0(substr(barcode, 1, 10), "T", substr(barcode, 11, 26)),
           paste0(substr(barcode, 1, 10), substr(barcode, 12, 26)))
  sequences <- paste0(prefix, bcs, suffix)
  records <- data.frame(
    name = c("normal", "insertion", "deletion", "reverse", "soft", "secondary", "supplementary", "partial", "hard"),
    flag = c(0, 0, 0, 16, 0, 256, 2048, 0, 0),
    pos = c(rep(1, 7), 21, 1),
    cigar = c("50M", "22M1I28M", "22M1D27M", "50M", "3S50M2S", "50M", "50M", "20S30M", "5H50M"),
    seq = c(sequences, sequences[1], paste0("TTT", sequences[1], "GG"), rep(sequences[1], 4)))
  lines <- vapply(seq_len(nrow(records)), function(i) {
    x <- records[i, ]
    paste(x$name, x$flag, "cassette", x$pos, 60, x$cigar, "*", 0, 0, x$seq,
          paste(rep("I", nchar(x$seq)), collapse = ""), sep = "\t")
  }, character(1))
  sam <- file.path(folder, "input.sam")
  bam <- file.path(folder, "input.bam")
  writeLines(c("@HD\tVN:1.6", "@SQ\tSN:cassette\tLN:50", lines), sam)
  stopifnot(system2(samtools, c("sort", "-o", shQuote(bam), shQuote(sam))) == 0L)
  stopifnot(system2(samtools, c("index", shQuote(bam))) == 0L)
  list(folder = folder, bam = bam, bcs = bcs)
}

test_that("flank extraction retains indels and excludes unsuitable alignments across chunks", {
  f <- flanked_bam_fixture()
  on.exit(unlink(f$folder, recursive = TRUE))
  out <- barbac_xtr(f$bam, "cassette", 13, 38, file.path(f$folder, "reads.csv"),
    flank_pattern = "GGTACC([ACGT]{24,28})ATAACT", include_read_ids = TRUE,
    yield_size = 1L, verbose = FALSE)
  x <- read.csv(out)
  expect_setequal(x$read_id, c("normal", "insertion", "deletion", "reverse", "soft"))
  expect_equal(x$barcode[match(c("normal", "insertion", "deletion"), x$read_id)], f$bcs)
  expect_equal(x$barcode_length[match(c("normal", "insertion", "deletion"), x$read_id)], c(26L, 27L, 25L))
  expect_equal(attr(out, "extraction_stats")$matched_alignments, 5)
  counted <- barbac_xtr(f$bam, "cassette", 13, 38, file.path(f$folder, "counts.csv"),
    flank_pattern = "GGTACC([ACGT]{24,28})ATAACT", yield_size = 2L, verbose = FALSE)
  counts <- read.csv(counted)
  expect_equal(counts$counts[match(f$bcs, counts$barcode)], c(3L, 1L, 1L))
})

test_that("query orientation and windows are explicit, and old extraction stays fixed width", {
  f <- flanked_bam_fixture()
  on.exit(unlink(f$folder, recursive = TRUE))
  out <- barbac_xtr(f$bam, "cassette", 13, 38, file.path(f$folder, "reverse.csv"),
    flank_pattern = "AGTTAT([ACGT]{24,28})GGTACC", reverse_complement = TRUE,
    include_read_ids = TRUE, read_window = c(1, 60), verbose = FALSE)
  x <- read.csv(out)
  expect_equal(x$barcode[x$read_id == "normal"], as.character(Biostrings::reverseComplement(Biostrings::DNAString(f$bcs[1]))))
  legacy <- barbac_xtr(f$bam, "cassette", 13, 38, file.path(f$folder, "legacy.csv"), verbose = FALSE)
  expect_true(all(read.csv(legacy)$barcode_length == 26L))
  expect_error(barbac_xtr(f$bam, include_read_ids = TRUE), "require flank_pattern")
  expect_error(barbac_xtr(f$bam, "cassette", 13, 38,
    flank_pattern = "(ACGT)", read_window = c(0, 26)), "read_window")
  expect_error(barbac_xtr(f$bam, "wrong", 13, 38, flank_pattern = "(ACGT)"), "Reference not found")
})
