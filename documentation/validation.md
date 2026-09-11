# Verification and benchmark evidence

## Publication candidate · 11 September 2026

The [release receipt](../benchmark/validation/release_2026-09-11/validation.json)
checks barbac 0.2.0 with the native **v14** clustering engine. The engine's R
clustering wrapper and C++ sources are unchanged from the validated development
revision `b7ce05c`; the package pipeline adds R1-only and mixed-input processing.

| Check | Result |
|---|---|
| Package tests | 44 cases, 526 passing assertions; no failures, errors, skips or test warnings |
| Studio integration | 8 cases, 57 passing assertions; no failures, errors or skips |
| R1-only and mixed inputs | Real FastQC, PEAR, minimap2 and samtools; all expected indel barcodes and counts recovered |
| Known-truth pipeline fixture | 144 reads per sample; three expected clusters; counts conserved |
| Studio extraction and plots | Single/paired extraction, native membership and plot parity, complete sample counts |
| Quarto exports | Self-contained HTML report containing the actual settings and results |

Run the current release verification from the repository root:

```sh
Rscript tools/verify_release.R
```

This requires the package and Studio dependencies, the four CLI tools above,
Quarto, `pkgload`, `testthat`, `jsonlite` and `digest`. It writes a separate receipt
under `benchmark/validation/release_2026-09-11/`. Timings are not speed benchmarks.
The local installed dependencies emit package-version and stack-imbalance
warnings when loaded; the tests themselves completed as recorded above.

## Archived paired-read verification · 10 September 2026

The [reproducible integration check](../benchmark/validation/verify_pipeline_extraction.R)
uses real FastQC, PEAR, minimap2, samtools and MultiQC executables, followed by
`barbac_xtr()`, clustering and R plots. It generates its own reads with known
barcode identities and counts; it does not download a sequencing dataset.
The dated [verification receipt](../benchmark/validation/pipeline_extraction/validation.json)
records the checked R/C++ source hashes, tool versions, test results and output hashes.
This archived run passed on **10 September 2026**, using **R 4.5.1 on macOS arm64**.

| Check | Verified result |
|---|---|
| FASTQ input and merging | 4 FASTQs, 2 samples, 290 pairs; 290 merged reads |
| Mapping and QC | 285 mapped reads and 5 deliberately unrelated unmapped reads; FastQC totals reconcile and MultiQC renders |
| Fixed-coordinate extraction | All 140 reads in the substitution-free sample have the exact expected barcode and count |
| Flank extraction | All 285 mapped reads have the expected barcode and read ID; observed 25-, 26- and 27-base sequences retain the designed indels |
| LV clustering and time series | 3 expected centroids; all 285 extracted reads conserved; all 6 lineage/timepoint counts match the known truth |
| R diagnostics and plotting | Extraction panels, mapping plot and `palette = "alger"` area plot render; frequencies sum to 1 at each timepoint |
| Package tests | 40 test cases, 488 passing assertions, 0 failures/errors/skips/test warnings |

Inspect the [barcode counts](../benchmark/validation/pipeline_extraction/barcode_counts.csv),
[lineage counts](../benchmark/validation/pipeline_extraction/lineage_counts.csv),
[extraction diagnostics](../benchmark/validation/pipeline_extraction/extraction_diagnostics.pdf)
and [time-series plot](../benchmark/validation/pipeline_extraction/barcode_time_series.pdf).
This is a small integration check on the recorded local environment. The
larger studies and accuracy benchmarks below provide separate evidence for
performance on experimental and simulated data. Docker and cluster execution
are not established by the local integration check.

## Larger applications and accuracy benchmarks

The [Chen 2023 application](../benchmark/time_series_chen2023/results/README.md)
processed 16,511,755 input pairs and retained 16,051,344 molecules after its
study-specific extraction and UMI filtering. The
[verification record](../benchmark/time_series_chen2023/results/validation.json)
checks input files, counts and normalization. Its
[workflow and SLURM instructions](../benchmark/time_series_chen2023/README.md)
and [R/Quarto report instructions](../benchmark/time_series_chen2023/r_report/README.md)
are available.

The [Jasinska 2020 application](../benchmark/time_series_jasinska2020/README.md)
processed 312 FASTQs containing 364,267,911 unique input reads, including the
shared initial samples. The [processing verification](../benchmark/time_series_jasinska2020/results/processing_validation.json)
reconciles all 78 sample receipts and FastQC totals. Its
[biological report source](../benchmark/time_series_jasinska2020/r_report/report.qmd)
compares three chloramphenicol-treated populations with three controls.
These archived analyses retain their original source/build hashes and timings;
they are distinct from the current branch integration check. Agreement with
published counts measures agreement with that processing reference, not known-truth accuracy.

The [paper comparison](../benchmark/latest_four_conditions/README.md) compares both
barbac modes with Shepherd, Starcode and Bartender on four simulated conditions
and the Milos/Johnson reference. These are **recorded v13 measurements**.
Native v14 preserves assignments in its recorded validation while improving LV
candidate search; no v14 timing is substituted into the v13 table.

Each cell below is **centroid F1 (%) / workflow time (seconds)**. F1 measures
recovery of true centroid identities: precision = TP/(TP+FP), recall = TP/(TP+FN),
and F1 is their harmonic mean. FP counts extra centroids; FN counts missed true
centroids. F1 is not the percentage of individual reads assigned correctly.
Accuracy averages three seeds in each smaller simulation; each workflow timing
is one serial observation. Milos uses one fixed input. R-S/R-I are random designs
with substitutions only/with indels; A-S/A-I are the corresponding anchored designs.

| Method | R-S | R-I | A-S | A-I | Milos |
|---|---:|---:|---:|---:|---:|
| barbac Hamming | 99.538 / 4.98 | 26.326 / 4.78 | 99.249 / 4.05 | 19.410 / 5.77 | 99.72147 / 20.81 |
| barbac LV + Poisson | 99.538 / 4.14 | 97.736 / 5.96 | 99.244 / 4.63 | 94.801 / 8.50 | 99.72246 / 43.40 |
| Shepherd | 99.502 / 3.47 | 79.648 / 3.74 | 99.162 / 117.94 | 65.188 / 105.08 | 99.72196 / 110.95 |
| Starcode sphere | 99.475 / 2.94 | 97.456 / 17.52 | 97.039 / 4.44 | 91.962 / 26.24 | 99.43682 / 155.30 |
| Starcode MP | 97.767 / 3.12 | 91.044 / 18.84 | 96.251 / 5.32 | 84.837 / 26.51 | 99.45408 / 150.22 |
| Bartender | 99.547 / 1.57 | 27.916 / 3.05 | 98.832 / 2.17 | 20.153 / 7.90 | 99.39120 / 36.52 |

A barbac mode has the highest reported F1 in three of four simulated conditions:
LV leads both indel conditions and Hamming leads A-S. LV also leads on Milos;
Bartender narrowly leads R-S. Speed rankings vary by dataset.
The [complete tables](../benchmark/latest_four_conditions/paper_table.md) include
FN, FP and assignment metrics, with the
[dataset schematic](../manuscript/media/benchmark_datasets.png) and recorded settings.
Earlier experiments remain in [the reference comparison](../benchmark/reference_comparison/README.md),
[the indel experiments](../benchmark/indel_experiment/) and
[the LV optimization notes](../benchmark/lv_optimization/README.md).

## Re-run verification

The repository contains code, reference sequences, compact summaries and small
verification artifacts. Raw sequencing files, full barcode-count matrices,
rendered HTML reports and their large plot assets are generated locally and
excluded from version control. Follow each study's scripts to recreate them;
paths inside report sources refer to those local outputs.

To rerun the earlier paired-read fixture, with the existing `barbac_env` available:

```bash
Rscript benchmark/validation/verify_pipeline_extraction.R
```

This older script replaces its own receipt when rerun; use an isolated checkout
to retain the archived record. The script needs `pkgload`, `testthat`, `jsonlite` and `digest` in addition to
barbac's dependencies. It runs the real paired-read pipeline, both extraction
modes, read-ID checks, clustering, plotting and the full package test suite.
It writes the receipt and small CSV/PDF artifacts under
`benchmark/validation/pipeline_extraction/`; temporary synthetic FASTQs and BAMs
are removed. A failed assertion exits with an error and writes a failed receipt.

For the package tests alone:

```r
devtools::load_all(".")
testthat::test_dir("tests/testthat", reporter = "summary")
```

