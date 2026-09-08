# Report design and validation notes

Audience: technical. Delivery: one portable HTML report, with its canonical
JSON payload and executable benchmark sources. The repository README documents
reproduction; it is not a second report implementation.

The technical-report specification is mapped as follows: title; result summary;
metric definitions moved ahead of evidence to prevent confusion about F1 versus
read accuracy; barcode recovery; runtime chart and component table; read
assignment; data-quality limitations; fixed methods; next steps and further
questions. No required section was omitted.

Chart contract: compare operational runtime across ten method configurations
on one fixed 100k-truth simulation. One horizontal bar per configuration,
sorted by seconds, zero-based scale; shorter means faster. Use a single palette
root because method identity is already on the categorical axis. Adjacent
prose states timing boundaries and the single-run limitation. The full reviewed
row retains FN, FP, TP, F1, read counts, conservation, and core/process timing,
so the source can support other comparisons. There is only one chart; exact
accuracy differences are better shown as counts in tables than as nearly
identical 99.99% bars. Time trends and distributions across replicates cannot
be drawn because there is only one supplied dataset and one timing observation
per configuration. Tables use explicit sorting appropriate to the question.

Data-quality audit: the companion notebook preserves the read-only profile.
The zero-read truths are an expected simulation condition with a large effect
on historical FN interpretation. Four inconsistent parent totals are a small,
high-confidence source issue; their cause is not established here. Label-based
accuracy is qualified and also recomputed excluding all reads belonging to
those four parents. This check does not validate the original simulation's
entire label-generation process. Temporal freshness/drift is not applicable.

Method validation: source/input multiplicity reconciliation, positive observed
counts, nonnegative truth counts, unique sequence keys, valid parent foreign
keys, unique predicted sequence assignments, and exact member-to-centroid count
reconciliation. The metric unit tests also verify that missing reads stay in
the denominator and invalid/duplicate mappings fail. All raw method outputs
are hashed and retained in the ignored work directory.

Timing scope: benchmark subprocesses ran sequentially. Lightweight source
inspection, report preparation, and small metric checks occurred during some
runs. This is an interactive development-machine experiment, not an exclusive
hardware timing study. Do not report uncertainty bounds or treat close timing
differences as established speed rankings.

Final report QA: canonical artifact validation and exact portable payload/structural verification passed. Chromium headless-shell was unavailable, so browser layout/source-dialog QA was not performed. The report retains semantic chart-data and table fallbacks. SQL recomputation of report metrics agreed with the Python scorer to 1e-14. Four independent metric tests pass.
