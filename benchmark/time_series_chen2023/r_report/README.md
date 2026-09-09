# Interactive report using barbac's R functions

Open **[report.html](report.html)** in a browser. It is a self-contained Quarto
report: all scripts, styles, figures and barcode count tables are embedded. No
R server, Python process, Shiny deployment or internet connection is needed to
read it. The complete counts make the file approximately 18 MB.

The report calls the package functions directly:

- `summarise_bam_stats()` recounts each of the eight completed BAMs;
  `plot_bam_stats()` supplies the mapping plot.
- `barbac_xtr.stats()` supplies all 16 sample/component extraction panels,
  including its length, abundance, entropy and summary-table output.
- `cluster_stats()` summarizes pooled component clusters and each sample's
  positive paired-lineage counts. Explicit zeros are retained in count tables
  and excluded only when counting observed clusters.
- `barbac_ts_area(interactive = "ggiraph")` supplies the four LV/Hamming ×
  replicate composition plots. The four largest pooled LV pairs plus Other
  preserve the full molecule denominator. The same IDs/colors appear in all
  four panels; `fill_missing = "zero"` preserves zeros.
- `theme_barbac()` supplies the shared plot theme. R's Plotly adapter adds
  interaction to mapping and comparison plots; DT supplies search and CSV export.

The numerical comparisons are recalculated in **R**, checked against all 48
existing sample/method results, and reported with their original caveats. The
report retains every barcode pair and all eight count columns for each barbac
mode. Its search and CSV buttons work entirely in the browser. CSV export uses
all filtered rows, not just the current table page.

## Rebuild with the existing environment

Run from the project checkout containing the completed mapped Chen workflow:

```bash
Rscript benchmark/time_series_chen2023/r_report/build_report.R \
  /absolute/path/to/chen_work
```

With no argument, the work directory defaults to
`benchmark/time_series_chen2023/generated/time_series`. The script loads the
current project with `devtools::load_all()` and reuses `barbac_env` for samtools.
It requires the existing `results/` comparison files, complete BAMs, extracted
component counts, and first-repeat centroid outputs. It writes only report
outputs and a separate cache under `generated/r_report`; it does not recluster,
re-extract, or change the previous result files.

The R reporting dependencies are `rmarkdown`, `knitr`, `htmltools`, `htmlwidgets`,
`ggiraph`, `plotly`, `DT`, `ragg`, `digest`, `jsonlite`, `base64enc`, and `devtools`,
plus the package dependencies. The Quarto CLI is required (local version 1.4.552). These were already installed
for the local run; no new environment or global package installation was needed.

## Files and validation

- `report.qmd`: Quarto narrative, package-function calls, and HTML widgets.
- `report_helpers.R`: checked data adapters and statistics preparation in R.
- `tables/`: numeric mapping, extraction, cluster and agreement summaries.
- `figures/`: the package's extraction panels and static composition plots.
- `validation.json`: input/source hashes and numerical check receipt.
- `qa/browser_validation.json`: browser control, offline, download and layout checks.
- `session_info.txt`: R and package versions used to generate the report.
- `plan.md`: function inventory, metric definitions and chart contracts.

Run the small adapter checks without sequencing files:

```bash
Rscript benchmark/time_series_chen2023/r_report/test_report.R
```

The package regression test is `tests/testthat/test-extraction-stats.R`. It
demonstrates a pre-existing `barbac_xtr.stats()` defect: the data column named
`barcode_length` masked the requested bounds inside `mutate()`. Explicit
`.env$barcode_length` now makes the summary table use the requested range,
independent of row order. The function signature, histograms, extracted counts
and clustering behavior are unchanged.

The original analysis is preserved at commit `bf171db`; its historical source
hashes are intentionally not rewritten to claim that this later QC fix was
used then. This report verifies those original result artifacts and records its
own R source hashes separately.

Optional browser QA uses the installed Chrome browser plus the R packages
`websocket`, `curl`, `later` and `jsonlite`. Start a temporary headless Chrome
with remote debugging on port 9223, then run `browser_checks.R`. It blocks all
HTTP(S) requests while loading the local HTML, checks every table and widget,
exercises sample expansion and composition tabs, and verifies an actual filtered
CSV download. Browser QA code evaluates JavaScript only to inspect and exercise
the browser UI; the report calculations and rendering remain R code.

Publication agreement is not ground-truth accuracy; 83.09% is the median LV
molecule coverage on exact published pair IDs. Clustering timing excludes the
shared FASTQ/BAM/extraction workflow and report generation. No fitness model is
fitted. See the report's definitions and limitations before reusing its numbers.
