# Reference dataset comparison

This experiment reruns the existing 100,000-barcode simulation that the project
previously compared with Shepherd, Starcode, and Bartender (the dataset referred
to as “Milos” in the discussion). It does not regenerate or modify the source data.

Inspect the [results CSV](results.csv). The generated `report.html` is a local
output; rendered reports and full data are excluded from version control.
Hamming with support ordering takes 18.4s for the measured workflow, versus
111.0s for Shepherd, with nearly identical recovery. Shepherd narrowly leads
F1 and supplied-parent agreement. Current LV with sequence ordering regresses
from 37.7s to 162.8s inside barbac compared with the previous implementation.

## Reproduce

Install the current R package into an isolated library, then run:

```sh
R CMD INSTALL --library=/path/to/current-library .
python3 benchmark/reference_comparison/run_comparison.py \
  --source /path/to/barbac-benchmark \
  --tools /path/to/barbac-benchmark/tools \
  --library /path/to/current-library \
  --baseline-library /path/to/previous-library
```

The Python environment needs pandas, numpy, and rapidfuzz. The previous library
is revision `0cbe94c`; the current clustering source is revision `dd6aad3`.
Without a previous library, select only current and external methods using
`--methods`. Tools and installed package hashes are recorded in `provenance.json`.

Input files: `barbac_benchmark_input.csv` (`barcode,counts`),
`simulated_reads.csv` (`BC,Count,true_BC`), and `true_counts.csv`
(`BC,True Count`). Inputs are staged in descending observed count and ascending
sequence order. No parent labels or true counts enter any method. All methods
receive the same sequences and multiplicities, including off-length sequences.
Bartender requires expansion to one row per read with unique read identifiers.

Large outputs default to the ignored
`benchmark/four_condition_comparison/generated/reference_2026-09-08/` directory.
Pass `--work` and `--report` for a new experiment. `--resume` validates saved
output hashes and reuses completed runs. Use a fresh work directory after
changing any library, tool, parameters, or input. Runs are sequential, with
one thread requested where supported and BLAS/OpenMP threads limited to one.

## What the scores mean

- FN: true barcode sequences absent from the output centroid set.
- FP: distinct output centroid sequences absent from the truth set.
- F1: `2 TP / (2 TP + FN + FP)`, with `TP = 100000 - FN`.
- Positive-truth F1/FN: the same calculation restricted to the 99,591 true
  barcodes with positive counts. Zero-read truths cannot be recovered from data.
- Observed FN: missing true barcodes whose exact sequences occur in the input.
- WS: FP centroids within Levenshtein distance three of any true barcode;
  this is a subset of FP and is not added again to F1.
- Read-assignment accuracy: sum of input counts whose predicted centroid exactly
  equals the supplied `true_BC`, divided by **all** 24,996,128 input reads.
  Wrongly assigned and unassigned reads both reduce the score.
- Consistent-parent accuracy: a sensitivity check excluding all source rows
  belonging to the four parent barcodes with inconsistent count totals.

The runner asserts that every reported member maps to only one centroid and
that summed member counts exactly reproduce reported centroid counts. It
checks missing reads explicitly. Bartender can generate consensus centroids
that never appear in the input; centroid scoring permits this.

## Timing

`core_seconds` is barbac's `super_cluster2` call, including input reading,
ordering, clustering, and result construction. `process_seconds` is an external
clock around the tool process, including startup and its native output.
`pipeline_seconds` additionally includes required method-specific preparation
and conversion to the common centroid/member format. It excludes shared input
staging, correctness evaluation, and hashing. For an operational comparison,
use pipeline time; compare core time only among barbac configurations.

One run per configuration is exploratory evidence, not a timing confidence
interval. Hardware, OS, commands, versions, and hashes are saved. Historical
numbers may differ because input tie ordering, versions, output requirements,
and timing boundaries differ. No algorithm parameters were tuned on this run.
There was also an initial Shepherd execution before fixing the wrapper's
output-path handling. Its process timing was not retained; its centroid and
membership outputs are identical to the reported execution. `PYTHONHASHSEED`
was not fixed, so no general claim of deterministic Shepherd tie handling is
made. See `shepherd_sensitivity.json`.

## Validation

```sh
python3 -m unittest discover -s benchmark/reference_comparison -p 'test_*.py'
```

The metric tests independently verify the treatment of absent truth, wrong and
unassigned reads, duplicated assignments, and count reconciliation.
`data_quality.ipynb` reproduces the source checks and displays their saved result.
`data_quality.json` preserves the per-parent discrepancies rather than editing
labels. Read-assignment accuracy measures agreement with these supplied labels;
it is not a claim of independently validated read-level ground truth.

`build_report.py` independently recomputes the report's accuracy metrics through
`report_metrics.sql` and asserts agreement with the Python evaluator. It writes
the portable report's canonical `artifact.json`; the Data Analytics plugin's
`skills/build-report/scripts/deliver_portable_artifact.mjs` packages it as HTML.
HTML payload and structure validation passed. Browser layout and source-dialog
interaction were not checked because the report builder found no installed
Chromium headless-shell. The HTML includes a readable semantic fallback.
