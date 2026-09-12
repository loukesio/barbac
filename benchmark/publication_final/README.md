# Publication benchmark: preserved v14 plus faster startup

This branch combines released v14 clustering with the independently verified
lazy BAM namespace loading change (`e91866a`). The clustering C++ and R wrapper
match main byte-for-byte. The tested installed library is reused after hash
verification. No calibrated, learned, or paired-indel experimental merge rule
is included. The R1 pipeline, app, extraction, palettes and plotting code remain
present in this release candidate.

## Frozen final evaluation

The final protocol registers **60 independent libraries per design**, random
N20 and anchored 26-base barcodes carrying the same 20 variable bases, plus the
unchanged published Milo reference. Each new library has 10,000 true identities,
one million expected reads, the existing abundance mixture, 0.004 substitutions
per base and the unscaled published homopolymer error table. Each method sees
identical input within a library.

Six configurations are measured: v14 Hamming, v14 LV with its existing Poisson
option, Shepherd, Starcode sphere, Starcode message passing, and Bartender.
Distances, ratios, ordering, nominal lengths and worker timing boundaries
retain their established settings. Hamming includes its existing rare-indel
rescue; it is not a strictly substitution-only comparator.

All 726 registered timed cells run serially with one thread per tool. This
includes one final-session Milo measurement per tool: the previous timings
predate the startup improvement and remain preserved. The explicit final
campaign is measured once; completed cells are never rerun for candidate tuning.
Inputs, source, binaries, commands, mappings, counts and receipts are retained.

The five previously reserved seeds are included, plus 55 newly registered
seeds. No final seed is used in the variance pilot or in model selection.
Failed cells remain failures. No seeds, rates, tools or thresholds are changed
after examining final results. The [failure-reporting addendum](FAILURE_REPORTING.md),
introduced after Shepherd's first automatic error-estimation failure, explains
how unavailable comparisons are represented. It preserves all registered cells
and limits statistical tests to contrasts with every required pair present.

## Replication and statistical claims

The planning pilot used four already-seen development libraries per design.
Existing results were reused; only 24 missing external-tool pilot cells were
run. The predeclared rule targets detection of a 0.01 percentage-point F1
difference, uses paired variability, requires at least 30 libraries, and caps
the campaign at 60 per design. The cap is binding: the pilot does not guarantee
90% power for every difference that small. Its four-library variance estimates
are uncertain. `sample_size.json` reports the uncapped requirements, detectable
effects and approximate power; the final sample size will not be extended in
response to observed results.

The eight primary contrasts compare LV with each external competitor in each
design, using **whole independent libraries as replicates**. Exact-identity F1
includes zero-read true identities, matching the original primary metric.
Simultaneous one-sided lower confidence bounds use a paired t approximation
with Bonferroni alpha/8. One-sided Holm-adjusted p-values are also reported.
A paired library bootstrap supplies a robustness check; a superiority label
requires both lower bounds to exceed zero. Parametric and bootstrap uncertainty
depend on their assumptions and finite replication; these are not universal
performance guarantees.

Positive-read truth F1, FN/FP, read assignments, abundance error and runtime
remain separate outcomes. Runtime intervals are secondary descriptive paired
log-time comparisons. Milo was used during development and is a single fixed
reference; it is excluded from independent-library significance claims.

## Explicit simulation boundary

The published calibration table covers repeat lengths 5–13. Our recursive
simulator updates coordinates and repeat lengths after each event. The original
archived notebook recursively copied its original `runs` descriptors; our
implementation is an adaptation, not byte-for-byte recreation of that notebook.

When a simulated read state requires an unmeasured repeat-length rate, it
becomes terminal for **further indel recursion**. All of those reads remain in
the dataset, undergo substitutions, retain their origins and are recorded in
`boundary_labels.csv`. Per-library counts show how often the boundary is used.
This finite-support modeling assumption does not imply that biological errors
stop at length 14. It avoids inventing unsupported rates or selecting libraries
that happen not to encounter the boundary.

Nine simulation/inference tests pass, including read conservation at the
boundary and paired statistical checks. Full regeneration of both original
mixed development inputs reproduces all four original input/truth files
byte-for-byte with zero boundary reads. Existing package validation remains
526 passing assertions on the same hash-verified speed library.

The previous stronger-rate pilot and its seven unsupported cells remain in
the paired-indel experiment record. Those 3x and 10x simulations are not silently
recast as empirical primary data. Genuine-neighbour controls and unsuccessful
merge candidates remain part of the supplementary evidence.

## Reproduction

Paths in `common.py` identify the preserved local references, tools and tested
library. Large generated artifacts are ignored by git. Compact protocols and
receipts are tracked.

```sh
python3 -m unittest discover -s benchmark/publication_final -p 'test_*.py' -v
python3 benchmark/publication_final/pilot.py
python3 benchmark/publication_final/run.py freeze
python3 benchmark/publication_final/run.py generate
python3 benchmark/publication_final/resume.py
python3 benchmark/publication_final/analyze_available.py
python3 benchmark/publication_final/report.py
```

The exact execution fingerprint is created after committing the protocol and
code and before any final generation. Keep frozen inputs and measurements;
do not amend them to fit later implementation changes.

The original launcher omitted a standard-library `json` import and stopped
after generating all inputs, before the first tool call. `resume.py` provides
that binding while preserving the original frozen script. Its separate
`execution_adapter.json` fingerprint is registered before measurement. No
generator, model, input, tool command, scoring or analysis changes accompany
this operational repair.

The original `analyze.py` is retained and refuses a campaign containing tool
failures. The separately recorded `analyze_available.py` reporting driver keeps
that frozen function's complete-pair calculations, fixed seeds and correction
family. Four additional checks verify failure propagation, unavailable means,
unchanged complete contrasts and preservation of the correction family.
`reporting_addendum.json` records when this addition was made and its hashes.
It is explicitly a post-start reporting decision, not retrospectively called
predeclared.

The main publication scope was chosen after development inspection; this
selection is disclosed. Background: [Johnson et al. barcode study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10276077/),
[archived analysis](https://zenodo.org/records/7411747), and
[Morris et al. simulation-study guidance](https://arxiv.org/abs/1712.03198).
