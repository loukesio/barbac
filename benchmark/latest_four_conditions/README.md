# Latest barbac comparison for the paper

The publication comparison uses **distance three**, five datasets, and six
configurations: barbac Hamming, barbac LV with the optional Poisson indel model,
Shepherd, Starcode sphere, Starcode default message passing, and Bartender.
Both barbac modes use native v13 and support ordering.

**Barbac has the highest centroid F1 in three of the four simulated conditions
and on the Milos reference dataset.** Bartender narrowly leads the random
substitution-only condition. Hamming leads anchored substitution-only data;
LV leads both indel conditions and Milos. This is a descriptive comparison of
the specified configurations, not evidence of universal superiority.

The four smaller conditions contain 10,000 true barcodes and one million reads
per seed, with results averaged over seeds 42, 43 and 44. The Milos dataset is
the Johnson et al. (2023) reference simulation used in the earlier benchmark:
100,000 listed true barcodes and 24,996,128 reads. FN and FP are counts; F1 is
exact centroid-recovery F1, not read-assignment accuracy.

Timings in the paper table are single serial workflow observations, using
seed 42 for each smaller condition. They include startup, required format
conversion, clustering, and both centroid/member exports. Earlier overlapping
peer timing measurements are excluded. The table does not establish timing
confidence intervals or statistical significance for small accuracy gaps.

## Paper deliverables

- [Table with caption](paper_table.md), [CSV](paper_table.csv), and
  [LaTeX](paper_table.tex).
- [Replacement results section](paper_results_section.md).
- Updated manuscript source: `manuscript/barbac_manuscript.md`, Sections
  3.3–3.5 and the comparison-dependent discussion.
- Generated editable Word file: `manuscript/barbac_manuscript_updated.docx`.
  This file remains ignored, following the repository's manuscript convention.
- [Arithmetic and scope validation](paper_validation.json) and
  [Word document validation](manuscript_validation.json).

The table retains both Starcode modes and the poor Hamming/Bartender indel
results. The optional LV indel model is explicitly identified. The reference
truth includes 409 zero-read barcodes; four per-parent count inconsistencies
are documented in the manuscript and the earlier reference audit. The smaller
simulations did not retain read-parent labels, so their read-assignment
accuracy is unavailable.

## Evidence and reproduction

`paper_source_runs.csv` contains the 78 observations used for accuracy:
72 simulation observations and six reference observations. `paper_table.csv`
contains the 30 publication rows. It records the accuracy sample size,
seed-specific timing basis, mean FN/FP/F1, and F1 range across seeds.

The runners use the original inputs under
`benchmark/four_condition_comparison/generated/revision-accuracy/`, the
validated isolated v13 library recorded in `manifest.json`, and the external
tools recorded in `completion_manifest.json` and `peer_manifest.json`.
Large inputs and membership exports remain in the ignored generated directory.
All reused inputs and outputs are hash-verified. Fresh barbac memberships
reconcile exactly with input counts and centroid totals. Saved Shepherd and
Starcode sphere accuracy results for seeds 43/44 are reused; seed 42 was rerun
serially. Both Starcode modes and Bartender have matching tool provenance.

The existing measured campaign is reproduced locally with these scripts in
sequence; the first runner retains the full exploratory campaign for audit,
and the paper builder selects only the requested distance-three comparison:

```sh
python3 benchmark/latest_four_conditions/run_comparison.py --resume
python3 benchmark/latest_four_conditions/time_peers.py
python3 benchmark/latest_four_conditions/complete_paper_comparison.py
python3 benchmark/latest_four_conditions/build_paper_table.py
python3 benchmark/latest_four_conditions/render_manuscript.py
```

To rebuild just the paper and Word document from saved measurements, run the
last two commands. The Word exporter verifies every table cell and all seven
original figure images, preserves the base document's embedded fonts, and
applies repeatable table headers. No clustering code changes were made for
this comparison.
