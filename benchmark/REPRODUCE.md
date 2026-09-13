# Reproduce the publication comparison

The archived campaign scripts and freeze files retain their original absolute
paths as execution evidence. Use `reproduce_publication.py` for a fresh checkout
on another machine; it accepts local tool/library paths and writes only to a new
output directory. It reuses the frozen simulator and scoring functions. The
completed Shepherd configuration is distance 3, Bayes threshold −4 and e=0.004
for new simulations, with automatic e for Milo.

The [0.2.1 release assets](https://github.com/loukesio/barbac/releases/tag/v0.2.1)
include all 120 canonical simulated libraries as `barbac-0.2.1-simulation-test-data.tar.gz`,
with per-file SHA-256 hashes. Use these to rerun the exact inputs without depending
on simulator byte reproducibility. The separate Milo reference remains linked below.

Prerequisites: Python 3.11+ with NumPy, pandas and SciPy; R and the Barbac package
dependencies. Original numerical-library versions are NumPy 2.2.6, pandas 2.3.3 and SciPy
1.15.3, recorded in `publication_final/freeze.json`. Install the checked-out release in an isolated
library, or use the library prepared by `Rscript app/run.R --prepare-only`. The
portable entry point reads Studio's active-library record by default; `--library`
can explicitly select another installation.

Generate one registered input without running any competitor:

```sh
python3 benchmark/reproduce_publication.py --seed 2026091202 \
  --condition random_mixed --output benchmark/reproduction/generated/random-1202
```

Evaluate an existing canonical input with Barbac only:

```sh
python3 benchmark/reproduce_publication.py \
  --source benchmark/reproduction/generated/random-1202/inputs \
  --condition random_mixed --methods hamming lv \
  --output benchmark/reproduction/generated/barbac-random-1202
```

For all six configurations add
`--methods hamming lv shepherd starcode_sphere starcode_default bartender`
and `--tools /your/tool/directory`. The tool directory must contain
`Shepherd/shepherd_t0.py`, `starcode/starcode`, and
`bartender-1.1/bartender_single`, built from the versions/fingerprints in the
archived receipts. A failed method is retained as a failure and is not retried.
No competitor runs occur unless explicitly included in `--methods`.

The complete set of 60 seeds is in `publication_final/final_protocol.json`.
Repeat generation for `random_mixed` and `anchored_mixed` with the same seed
to reproduce paired designs. Exact byte reproduction also depends on the
recorded Python/numerical-library environment. The final seeds are evaluation
data and must not be used for further model tuning.

For Milo, obtain [Johnson's archived analysis and data](https://zenodo.org/records/7411747)
and supply `--condition milos --source /path/to/canonical/milos`. The canonical
files are `input.csv` (`barcode,counts`), `truth.csv` (`barcode,true_count`), and
`labels.csv` (`member,true_barcode,read_count`). Their contents/hashes and the
original conversion are recorded in `publication_final/run.py` and the frozen
dataset metadata. Preserve all 100,000 truth identities, including zero-read
identities, and all 24,996,128 read labels.

Saved-table reproduction needs no new clustering:

```sh
python3 benchmark/shepherd_completion/publication_table.py
```

The table renderer additionally uses python-docx, ReportLab, the documented
Times New Roman fonts and PDF export utilities. Its provenance checks reconcile
all values with the saved execution records. Independent reruns produce new
runtime observations; they never overwrite the recorded publication times.
