# Four-condition clustering comparison

The newer [exact-search experiment](EXACT_SEARCH.md) compares separate installed
revisions, both barbac distances, and an optional support-based tie rule over
multiple simulation seeds. The September 1 tables below are retained as history.

To run a fresh paired comparison after installing each revision in its own R
library:

```bash
python3 benchmark/four_condition_comparison/compare_revisions.py \
  --baseline-library /tmp/barbac-baseline \
  --candidate-library /tmp/barbac-candidate \
  --seeds 42 43 44 --repeats 3
```

The runner records individual timings and output hashes, verifies reused inputs,
and supports `--resume` only for the same configuration and installed builds.
Use `time_revision.R` for repeated timings in one loaded R session; keep those
separate from end-to-end process time.

This directory preserves the code, configuration, tool revisions, and compact
results for the 2026-09-01 comparison of barbac, Shepherd, Starcode, and
Bartender.

## Conditions

Every condition uses 10,000 true barcodes, 1,000,000 reads, lognormal
abundances, simulator seed 42, and truth-independent sequence ordering for
equal-count observations.

| design | errors | barcode |
|---|---|---|
| random | 0.5% substitutions | 20 random bases |
| random | 0.5% substitutions + 0.5% insertions + 0.5% deletions | 20 random bases before errors |
| anchored | 0.5% substitutions | `NNNNNNNNATGCNNNNNNNNATCGTTAA` |
| anchored | 0.5% substitutions + 0.5% insertions + 0.5% deletions | same template |

Rates are per base. Pure-substitution conditions run barbac in both Hamming and
Levenshtein modes. Indel conditions run barbac in Levenshtein mode only.

## Reproduce

From the repository root:

```bash
python3 benchmark/four_condition_comparison/run_benchmark.py
```

The script:

1. clean-builds the current barbac source into a temporary R library;
2. verifies the expected native build identifier;
3. deterministically regenerates all four simulated datasets;
4. runs every method and evaluates FN, FP, WS, Pearson correlation, and time;
5. writes raw outputs, `summary.csv`, and `versions.json` under `generated/`.

`generated/` is ignored because it contains large reproducible intermediates.
The compact dated snapshots are tracked:

- [`summary_2026-09-01.csv`](summary_2026-09-01.csv)
- [`versions_2026-09-01.json`](versions_2026-09-01.json)

The external tools are expected under
`~/Documents/Projects/Barcodes/barbac-benchmark/tools/`, matching
`benchmark/indel_experiment/run_experiment.py`. Python requires NumPy, pandas,
RapidFuzz, and SciPy; R requires the package dependencies in `DESCRIPTION`.
Runtime is machine-dependent and the complete run may take roughly 10–15
minutes. Bartender temporarily expands collapsed counts to one row per read.

## Snapshot results

Cells show `FN / FP / WS; algorithm seconds`. Lower is better for all four
values.

| condition | barbac | Shepherd | Starcode | Bartender |
|---|---:|---:|---:|---:|
| random substitutions | Hamming: 54 / 58 / 54; 0.65s<br>LV: 54 / 58 / 54; 1.19s | 54 / 58 / 54; 4.41s | 48 / 51 / 47; 4.56s | 48 / 52 / 48; 1.63s |
| random substitutions + low indels | LV: 170 / 350 / 171; 2.92s | 103 / 4,924 / 4,890; 5.56s | 158 / 338 / 159; 30.83s | 94 / 51,102 / 50,925; 4.33s |
| anchored substitutions | Hamming: 79 / 78 / 69; 4.02s<br>LV: 80 / 78 / 69; 12.92s | 79 / 78 / 69; 190.37s | 319 / 266 / 257; 6.80s | 105 / 108 / 100; 2.23s |
| anchored substitutions + low indels | LV: 266 / 957 / 271; 38.99s | 177 / 10,340 / 10,199; 163.73s | 476 / 1,180 / 494; 39.80s | 206 / 77,608 / 76,934; 10.65s |

FN is the number of true barcodes absent from the reported centroids. FP is the
number of reported centroids not in truth. WS is the subset of false-positive
centroids lying within edit distance three of truth, representing likely
wrongly split error sequences. Low FN is not useful by itself when FP and WS
explode, as happens for Shepherd and Bartender in the indel conditions.

For barbac, algorithm time excludes R interpreter/package startup. Peer tools
have negligible startup relative to their reported process runtime. The LV and
peer-tool rows were produced with native build
`barbac-2026-09-01-composition-prefilter-v10`. The two Hamming rows were rerun
with `barbac-2026-09-01-hamming-refinement-v11`; anchored Hamming's complete
centroid set and every cluster count were identical to Shepherd. The exact
source worktree was dirty, so the tracked changes must be committed together
with this artifact before the Git revision alone becomes sufficient provenance.
