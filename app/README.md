# barbac Studio

A Shiny workspace for extracting DNA barcodes, clustering them with the native
barbac engine, exploring lineage trajectories, and downloading a reproducible
analysis. The app uses the code on this branch; no alternative clustering
implementation is maintained inside the interface.

## Watch the workflow

[![Explore barbac Studio results and lineage plots](media/studio-preview.gif)](media/studio-walkthrough.mp4)

| Video | Workflow |
|---|---|
| [Barcode counts → complete analysis · 44 seconds](media/studio-walkthrough.mp4) | Upload a count table, run native LV clustering, inspect trajectories, palettes and memberships, download the results ZIP, and generate an HTML report |
| [FASTQ → extracted barcodes · 24 seconds](media/studio-fastq.mp4) | Upload paired FASTQs and a reference, apply extraction settings, download extracted counts immediately, and continue to clustering |

The [offline video viewer](media/watch.html) plays both clips in a browser with
no R session or internet connection required. Keep it alongside the videos and
poster images in `media/`. The recordings show actual local app interactions
using only synthetic fixtures, with captions and no playback speedup. They are
workflow demonstrations, not performance benchmarks. The animated preview is a
ten-second excerpt of the first video. Recording details are in
[recording.json](media/recording.json).

## Start locally

From the repository root:

```sh
Rscript app/run.R
```

Open **http://127.0.0.1:3838**. The launch script builds the current package with
normal release compiler flags into `app/.runtime/library`, keeping it separate
from your global R library and the archived benchmark installations. It rebuilds
when package source hashes change. The first launch needs the package's normal
dependencies and a working C++ compiler.

Studio additionally uses these R packages:

```r
install.packages(c("shiny", "bslib", "DT", "future", "promises", "jsonlite",
                   "digest", "ggiraph", "zip"))
```

Shiny >= 1.8.1 is required. Install the [Quarto CLI](https://quarto.org/docs/get-started/)
for HTML report exports. Raw-read extraction additionally requires minimap2 and
samtools; paired reads also require PEAR. Studio resolves those executables from
barbac's existing environment or PATH. Use `barbac::configure_environment()` and
`barbac::check_barbac_tools()` when setting up a new machine.

## Extracted barcode input

Upload one or more CSV/TSV files (optionally gzip-compressed). Column names are
case-insensitive. Required columns:

```csv
barcode,counts
ACGTACGTACGT,1500
ACGTACGTACGA,12
TTGCAATTGGCC,830
```

Optional columns are `sample`, `time` and `population`. Without `sample`, the
filename supplies the sample label. Without `population`, all samples are one
library. Without timepoints, clustering and cluster statistics work, and the app
explains why no trajectory can be drawn. A separate metadata CSV can provide
`sample,time,population`; its sample labels must match the uploaded tables.

Use a different population label for each independent barcode library or
replicate series. Every sample must have one population and one numeric timepoint
(or no timepoints for that whole population). Two samples cannot silently share
a population/timepoint. Exact duplicate sequences are summed within each sample.
Counts must be positive whole numbers. DNA is normalized to uppercase; A/C/G/T/N
are accepted. Invalid rows cause an actionable error rather than being dropped.

The included example is **synthetic**: 18 parent sequences, eight one-substitution
variants, eight timepoints, and 192,000 total barcode reads. It is an interface
demonstration, not a biological result or performance benchmark. Example tables,
sample metadata, and a small synthetic FASTQ bundle are downloadable in the app.

## Raw FASTQ input

The app supports single-end reads and overlapping paired-end reads, with one
barcode locus and one reference cassette per run. For multiple paired samples,
R1 and R2 files must be supplied in matching order. The app checks every FASTQ
record, paired read counts, and paired identifiers before processing.

1. Paired reads are merged with PEAR. Single-end reads map directly.
2. Minimap2 maps reads to the supplied cassette; secondary/supplementary
   alignments are removed before extraction. Samtools sorts and indexes the BAM.
3. `barbac_xtr()` extracts the barcode using explicit reference coordinates and
   either exact flanks in query space or fixed reference coordinates.
4. Extracted counts can be downloaded immediately or used for clustering.

Exact-flank extraction preserves observed lengths and accepts separate minimum
and maximum lengths. The flanks are literal DNA, not user-supplied executable
code or arbitrary regular expressions. They must match the reference-oriented
query. Fixed-coordinate extraction does not preserve observed indel lengths.

This route counts **merged molecules** for paired data and excludes unmerged
pairs. It does not implement UMI deduplication, FastQC/MultiQC reporting, or a
study's specialized filtering. Matching a published analysis requires its own
documented preprocessing. The synthetic download includes a reference, exact
settings, expected counts, paired reads and a separate single-end example.

## Scientific behavior

`super_cluster2()` runs on pooled counts within each population, then memberships
are mapped back to every original sample. The adapter checks membership coverage,
population totals and every sample total. Pooling is retrospective: later
timepoints can inform the cluster identity used at earlier timepoints.

Defaults match the package: LV, distance 3, merge ratio 20, error rate 0.005,
sequence ordering, and no Poisson indel exception. Neighbour-support ordering and
the experimental Poisson option are explicit choices. Hamming in Studio is
restricted to fixed-length A/C/G/T sequences up to 32 bases; use LV otherwise.

`cluster_stats()` supplies the actual numerical summaries. `barbac_ts_area()`
draws all inferred lineages with `min_total_count = 0` and `fill_missing = "zero"`.
Frequencies divide by all supplied barcode counts in each sample. A palette change
does not recluster data. Dense plots above 1,500 lineages use a static rendering;
grids above two million lineage/timepoint cells are left to the complete CSV/R
export instead of freezing the browser. No remainder category is introduced.

Each completed run saves centroids, full memberships, sample counts/frequencies,
statistics, an RDS result, parameters, the native build ID, source fingerprint,
and normalized-input hash. Raw extraction also records input/reference hashes
and processing counts. HTML reports are rendered from `report.qmd` in a background
worker and embed their resources. Clustering timings exclude upload, extraction,
plotting and report generation.

## Local runtime

Studio is kept on the local computer. The default address, `127.0.0.1`, accepts
connections from that computer only. No public app has been deployed. Use the
README videos for a walkthrough without running the analysis environment.

Uploads and job directories are private to a Shiny session. Closing the session
removes temporary files; a running worker cleans them up when it finishes. Download
results before closing the browser. No uploads or results are committed to Git.

Defaults can be configured with environment variables:

| Variable | Default | Meaning |
|---|---:|---|
| `BARBAC_STUDIO_PORT` | 3838 | HTTP port |
| `BARBAC_STUDIO_HOST` | 127.0.0.1 | Bind address |
| `BARBAC_STUDIO_WORKERS` | 2 | Shared background R workers |
| `BARBAC_STUDIO_UPLOAD_MB` | 512 | Upload request and combined raw-read limit |
| `BARBAC_STUDIO_MAX_ROWS` | 2000000 | Maximum barcode input rows |

The app rejects population totals above the current native engine's signed
32-bit count capacity. These are admission limits, not memory guarantees.
Large diverse libraries may take tens of minutes and need enough local CPU/RAM.

## Verification

```sh
Rscript app/tests/run_tests.R
Rscript app/scripts/browser_checks.R
```

The browser checks expect Studio at port 3838 and Chrome's DevTools endpoint at
port 9223. They exercise real uploads, background clustering, native interactive
plots, palette changes, ZIP/HTML downloads, input errors, paired extraction, and
desktop/mobile layouts. Local receipts and screenshots are written to `.qa/`.
The engine tests independently compare app memberships with the native API,
check count conservation and population separation, reject invalid inputs,
recover known truth from real single/paired FASTQs, and inspect the Quarto export.

To record the walkthroughs again, use `Rscript app/scripts/record_demo.R` with
the same local app and Chrome endpoints. Recording also requires `ffmpeg` on PATH.
