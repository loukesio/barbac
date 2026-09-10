# Jasinska 2020: E. coli lineage time series

This is a biological application of barbac LV, separate from the multi-method
benchmark. The selected study is [Jasinska et al. (2020)](https://doi.org/10.1038/s41559-020-1103-z),
*Chromosomal barcoding of E. coli populations reveals lineage diversity dynamics
at high resolution*. The report uses R, barbac diagnostics, Quarto, gt tables,
ggiraph and Plotly. Counts represent single-end sequencing reads, not cells or
UMI-deduplicated molecules.

## Selected experiment

Selection was fixed before examining barbac agreement: the **constant-concentration
experiment**, BioProject **PRJNA592529**, with **chloramphenicol 1 µg/mL** and
**no antibiotic**, three biological replicate populations per condition, and all
their deposited timepoints. See [samples.tsv](samples.tsv) and
[sample-to-condition mapping](sources/well_conditions.csv).

- Low CMP: wells A3, B3 and C3; passages 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 30.
- No drug: wells A1, B1 and C1; the same passages plus passage 1.
- Four sequencing subsamples are combined into each sample. These are not
  additional biological replicates.
- Three initial-population samples, each with four subsamples, are pooled into
  one common baseline. See [baseline_samples.tsv](baseline_samples.tsv).
- 75 longitudinal samples / 300 FASTQs / 350,934,958 reads / 8.384 GB compressed,
  plus 12 baseline FASTQs / 13,332,953 reads / 0.330 GB compressed.
- The source table converts each passage to approximately six generations.

The initial library, increasing-drug experiment and other treatments are not
part of this selected application. The completed Chen yeast analysis remains
available in its own folder.

## What can be compared with the publication

The author supplement [supplementary_tables.xlsx](sources/supplementary_tables.xlsx)
contains Supplementary Table 1c (per-sample input, quality, extraction and barcode
richness) and Table 4b (the top 20 barcode identities and mean/final frequencies
per population). The parser preserves the published 13- and 14-base identities;
15 bases is the nominal construct length, not an excuse to drop shorter references.

These tables permit exact input-count reconciliation, a richness comparison,
and final-frequency comparisons for all 120 selected population/barcode entries.
They do not contain the complete published timepoint-by-barcode count matrix.
Do not describe the application as reproducing every trajectory, all diversity
indices, or every biological conclusion. Author mean frequencies are retained;
initial-sample and averaging conventions must be resolved before treating their
differences as a definitive reproduction test. The paper states that custom
analysis scripts are available on request. No messages to authors have been sent.

## Reference and extraction

[cassette.fasta](reference/cassette.fasta) transcribes the **288-base construct
sequence printed in the paper's Methods**, with Ns at positions 11–25 (one-based,
inclusive). It is a barcode-cassette reference, not a whole E. coli reference genome.
The nominal barcode is 15 random bases and the downstream sequence begins
`TATCTCGGTAGTGGGATACGACGATACCGAAGACA`.

The single-end reads start near the variable region. Minimap2 aligns their
downstream constant sequence using the short-read preset and `-k 9 -w 5 -m 10
-s 10 -n 1 --secondary=no`. Secondary/supplementary records are excluded before
the sorted BAM is counted. Barcode bases remain in the query, including soft
clipping; no masked-reference base is substituted into a barcode.

`barbac_xtr()` selects alignments covering the downstream anchor at positions
32–40 and captures `^([ACGT]{10,20})TATCTCGGTAG` from the query. This accepts
observed 10–20-base sequences and preserves their variable lengths. The anchor
coordinates select the cassette; they are not the barcode's reference coordinates.
This flank rule is explicit and differs from the authors' custom alignment rule.
Extraction discrepancies are reported, not attributed entirely to clustering.

The quality filter matches the paper's minimum Phred 10 over each entire read.
No UMI deduplication is performed. FastQC is run on each deposited FASTQ; quality,
primary mapping and extraction counts reconcile to the input manifest. Generated
receipts record input MD5, reference/source hashes, output hashes and timings.

## Clustering and diversity

Barbac LV settings are fixed in advance: distance 3, support ordering, merge ratio
20, assumed error rate 0.005, Poisson indel model, design scoring off. Each population
is pooled over all its timepoints and the shared baseline; every member assignment
is mapped back to each sample and read totals must be conserved. This is retrospective
reconstruction. The baseline is reused for clustering but is counted once in unique
sequencing totals. No competitor executable is run by this application.

Area-plot frequencies divide each barcode count by all extracted barcode reads
at that passage. Every inferred barcode is drawn separately, with no abundance
cutoff or combined remainder. A shared plasma palette assigns colours by barcode
order across all six populations. Explicit zero counts stay zero. Read-accounting
tables retain all-input denominators. The author final-frequency denominator is unresolved:
Low CMP r1's 20 published final frequencies sum to 86.94%, exceeding the 80.24%
of reads identified as barcodes in Table 1c at passage 30. Both cannot use all
input reads. This could reflect different normalizations or a mismatch between
source tables; it is not resolved by the available custom-code description.
`frequency_audit.R` therefore reports both all-input and extracted-read
normalizations, with both absolute differences and the scale-invariant rank
correlation. No normalization is chosen solely for favorable agreement.

The six area panels are R-generated images, with one polygon per lineage.
`barbac_ts_area()` prepares all frequencies and colours; a vectorised grid draw
renders the same completed stack efficiently for roughly 200,000 bands per plot.
Every band's width is checked against the source count matrix. Individual-band
hover is unavailable at this scale; the report's other interactive charts and
comparison table remain available. `r_report/barcode_colours.csv.gz` records the
shared mapping; full barcode count matrices remain in `results/`.

Population diversity uses full counts, never top-N display categories. Exports
contain richness, exponential Shannon and inverse maximum frequency under both
total-read and extracted-read denominators. The main diversity panels normalize
among extracted barcode reads, making standard effective-number definitions
valid. Do not confuse these with within-sequence nucleotide-composition entropy
in the extraction panels, or with the paper's per-position library entropy.

## Reproduce with the existing barbac environment

Run from the repository root with the existing R dependencies and `barbac_env`.
The analysis code is R; minimap2, samtools and FastQC are external command-line
tools already used by barbac. The report additionally needs gt, readxl, xml2,
digest, jsonlite, viridisLite, plotly, DT, ragg and the Quarto CLI.

```bash
Rscript benchmark/time_series_jasinska2020/prepare_sources.R
Rscript benchmark/time_series_jasinska2020/download_reads.R
Rscript benchmark/time_series_jasinska2020/build_release.R
Rscript benchmark/time_series_jasinska2020/process_samples.R
Rscript benchmark/time_series_jasinska2020/validate_processing.R
Rscript benchmark/time_series_jasinska2020/cluster_populations.R 3
Rscript benchmark/time_series_jasinska2020/benchmark_release.R 3
Rscript benchmark/time_series_jasinska2020/validate_results.R
Rscript benchmark/time_series_jasinska2020/r_report/build_report.R
```

Source snapshots required by `prepare_sources.R` are tracked in `sources/`.
Downloads are checked against ENA byte counts and MD5 hashes. A complete analysis
needs roughly 9 GB for compressed FASTQs plus space for BAMs and temporary
quality-filtered FASTQs; keep at least 25 GB free. Processing uses two samples
concurrently, with two mapping threads per sample. Do not launch overlapping
processing commands into the same sample directories.

The HTML report is self-contained for its figures, summary tables and published
dominant-lineage explorer. Complete reconstructed count matrices are separate
compressed CSVs in `results/`; they are not all embedded in the HTML.

```bash
Rscript benchmark/time_series_jasinska2020/test_analysis.R
Rscript benchmark/time_series_jasinska2020/test_extraction.R
Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-extraction-stats.R")'
```

The report's A–C panels and numeric summary come from the optional
`barbac_xtr.stats(..., panel_labels=TRUE, return_details=TRUE)` interface. The
default still returns its original patchwork object. The report formats the
returned length summary as native HTML with `gt` and explains both denominators.

## Cluster execution

The [SLURM script](run_analysis.slurm) runs processing, clustering and reporting
sequentially within one job. Stage the downloaded FASTQs and load the site's R,
Quarto and existing barbac environment before submission. Add your site's account
and partition options to the submission command as required:

```bash
sbatch --export=ALL,BARBAC_PROJECT=/path/to/barbac \
  benchmark/time_series_jasinska2020/run_analysis.slurm
```

This script has not been submitted to a cluster. It requests eight CPUs, 24 GB
RAM and a 12-hour time limit; actual resource use depends on the cluster and
whether outputs are already cached.

`cluster_populations.R` runs independent populations in separate R processes
(up to six on this SLURM allocation), checks every exit status, then combines
the completed caches. It verifies the isolated release installation before starting workers.
Outside SLURM it defaults to two workers; an explicit worker count can be supplied,
for example `Rscript benchmark/time_series_jasinska2020/cluster_populations.R 6`.
Do not launch another clustering command for the same wells while it is running.

## Exact-search optimization required by this dataset

The short, dense barcode library exposed a large candidate-search cost in v13.
Build v14 changes only the wider LV search: after the complete distance-one
search, a child with count at least five cannot merge at distance two or greater
into a parent below the distance-two abundance requirement. With the default
ratio, that lower bound is 60 times the child count. Larger distances have
stronger guards; the distance-one Poisson exception has already been checked.
The bound is disabled when a usable design mask can override the guard.

This skips candidates that cannot change the assignment. It changes the debug
counts of evaluated blocked candidates; it does not change the scoring or merge
rules. The full test suite checks indexed/exhaustive equivalence, including
indels, abundance-floor boundaries, unsorted native inputs and design cases.
[search_validation.json](sources/search_validation.json) records identical
centroids, member assignments and counts for 12,000 real input sequences under
v13 and v14. That small check is an equivalence test, not evidence of a general
speedup. The incomplete v13 full-population run is not a full timing benchmark.

Timing exports distinguish elapsed time from active process CPU time.
`devtools::load_all()` defaults to a debug build (`-O0`); the original local
clustering times are archived in `results/clustering_times.csv`. The report uses
`release_clustering_times.csv`, measured with a separately installed optimized
package (`-O2` on this machine). `benchmark_release.R` verifies every membership
and centroid count against the original results and records the compiler command,
binary hash and source hashes in `release_validation.json`. It reuses timings
only when the original clustering cache records the identical release binary.

`build_release.R` installs under `generated/release/lib` without replacing the
user's installed barbac. Build once before clustering; `load_release.R` rejects
stale sources or binaries. A fresh SLURM run uses this installation from the
start. Historical benchmark tables still describe their original recorded
builds and have not been rewritten as v14 measurements.

## Validation of the report

`validate_results.R` independently reconstructs every exported count-matrix
cell from raw barcode tables and the recorded cluster memberships, checks all
75 sample input totals against the publication, and verifies the diversity
calculations and all 120 published final-frequency entries.
`validate_processing.R` also reconciles all 312 independent FastQC read totals
and exports the 78 source-hashed sample receipts.

`r_report/browser_checks.R` uses R to drive a temporary local Chrome instance
on port 9223. It checks offline rendering, native gt tables, the A–D panels,
three replicate tabs showing treatment and control side by side, their six
full-barcode figure receipts, barcode search,
CSV download and a narrow
mobile viewport. JavaScript expressions in that file inspect the browser UI;
all biological calculations and report preparation are in R.

## Optional LTC palette comparison

The [palette comparison](r_report/palette_comparison.html) shows every barcode
from treatment replicate 1 beside control replicate 1 using `alger`, `dora` and
`casa_natal`, plus a swatch sheet of all 32 palettes exposed by
[`ggvmap::vm_palettes()`](https://github.com/loukesio/ggvmap). Frequencies, barcode
order and polygon geometry are identical across choices. The palette mapping is
shared between treatment and control. The main report retains plasma while these
alternatives are reviewed. No palette is claimed to make hundreds of thousands
of adjacent bands visually distinguishable.

The optional comparison uses the installed `ggvmap` package (recorded version
0.3.0); barbac itself already accepts its colour vectors through `palette=`.
To rebuild the previews after building the main report assets:

```bash
Rscript benchmark/time_series_jasinska2020/r_report/palette_preview.R swatches
Rscript benchmark/time_series_jasinska2020/r_report/palette_preview.R A3
Rscript benchmark/time_series_jasinska2020/r_report/palette_preview.R A1
Rscript benchmark/time_series_jasinska2020/r_report/build_palette_comparison.R
Rscript benchmark/time_series_jasinska2020/r_report/browser_checks.R
```

Main-report figures use separate R processes for parallel rendering. Forked
rendering can fail during macOS font initialization; failed workers must stop
the build instead of producing a success receipt.

Checksum maps are written as JSON objects keyed by filenames. Readers reject
empty or unnamed maps, missing files and changed contents. The local
`metadata_schema_migration.json` receipt documents a one-time correction of
name-dropping JSON serialization, with all six original cache signatures and
timings preserved. `repair_metadata.R` is only for the saved local pre-migration
receipts; fresh analyses use the corrected schema directly.
