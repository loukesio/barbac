# Handoff: barbac clustering performance and benchmarking

**Complete Chen 2023 time series (9 September 2026):** all eight hBFA1/YPD
samples (generations 8/16/24/40, two biological replicates) have run locally
through FastQC → PEAR → minimap2 → BAM → indel-preserving `barbac_xtr()` →
original-order quality/UMI filtering → six-method clustering. Reused
`barbac_env`; no real SLURM jobs were submitted. The empirical 167-base cassette
reference has BC2 at 50–75 and BC1 at 110–135 (one-based; reverse-complement BC1).
Observed 24–28-base components and paired identities are retained.

See [results and figures](benchmark/time_series_chen2023/results/README.md),
[reproduction/SLURM instructions](benchmark/time_series_chen2023/README.md), and
[validation](benchmark/time_series_chen2023/results/validation.json).
All 16 complete FASTQs pass ENA MD5/size checks; FastQC reconciles 16,511,755
input pairs. Extraction retains 16,051,344 molecules (97.21%). Both pooled
components are clustered at distance 3; barbac uses support ordering, ratio 20,
error rate 0.005, and Poisson for LV only. Native clustering remains v13.
Median Spearman agreement with the 2,314 published pair IDs is 0.999968 for LV
(14.65 s combined-component workflow) and 0.999929 for Hamming (11.97 s).
Starcode has slightly higher reference agreement; Hamming is fastest by median.
All memberships match across three fresh serial repeats. Timing excludes shared
preprocessing. Shepherd automatic error estimation failed on BC1; both its
components use documented `-e 0.005`, set before publication comparison.

Median LV molecule coverage on exact published IDs is 83.09%; this is not an
accuracy percentage. Pooled unmatched mass is 19.70%, of which two abundant
pairs account for 69.47%; their BC2 sequences are at least seven edits from every
published BC2. The reason the author reference omits them is unresolved.
Publication counts are a differently processed reference, not true FP/FN labels.
Section 3.6, Table 3 and Figure 6 document this application; Table 2 retains the
previous five-dataset benchmark. The Word manuscript is regenerated locally.

The optional flank mode in `barbac_xtr()` adds query-space extraction without
changing default fixed-coordinate behavior. Synthetic indel/orientation tests,
13,657-alignment real-BAM parity, and ten Python extraction/submission tests
pass. Full `R CMD check --no-manual` has zero errors and one existing AppleClang
warning from R's `R_ext/Boolean.h`; tests and vignettes pass. Reference-construction
and direct-parser pilot reports are historical validation records. Cluster
account/partition/scratch configuration remains site-specific.

The five-dataset comparison now also has a compact F1/time summary and a
publication schematic in PNG/SVG/PDF, embedded as Figure 5 in the regenerated
Word manuscript. Numerical benchmark measurements are unchanged.

**Latest paper update:** See
[the five-dataset comparison](benchmark/latest_four_conditions/README.md).
The requested comparison uses distance three, both current barbac modes,
Shepherd, Starcode sphere and message passing, and Bartender. Table 2 in the
manuscript now contains all 30 dataset/configuration rows, with three-seed
accuracy means for the four simulations and fresh serial seed-42 workflow
timings. Milos uses the full validated reference dataset; Hamming was rerun
with v13 (20.81s). Barbac leads centroid F1 in three simulation conditions and
on Milos; Bartender narrowly leads random substitutions. The Word manuscript
is regenerated locally from the tracked Markdown. The user explicitly wants
this method comparison, not a distance-two-versus-three report.

**2026-09-08 update:** Work now continues on `feat/exact-search-clustering`.
**Latest: native v13** adds exact abundance-bound pruning and an opt-in
`indel_model = "poisson"` for repeated-base single indels. Read the
[LV optimization experiment](benchmark/lv_optimization/README.md): LV support
now takes 43.4s including exports, with FN 471 / FP 83 and 131 wrong read
assignments, versus the v12 values 192.0s, FN 471 / FP 85, and 231 wrong reads.
Search-only v13 preserves all reference member assignments. Across four
conditions and four seeds, the optional model adds no FN and removes six FP.
The homopolymer-rich stress case improves but remains weak. The indel model
stays off by default; no universal accuracy or speed claim is established.

The earlier [100,000-barcode reference comparison](benchmark/reference_comparison/README.md)
tests both current distance modes and ordering options against the previous
version, Shepherd, both Starcode modes, and Bartender. Current Hamming with
support ordering has FN 469 / FP 87 in 18.4s including startup and exports;
Shepherd has FN 470 / FP 85 in 111.0s including output conversion and narrowly
better F1 and supplied-label agreement. Current LV with sequence ordering has
the same FN/FP as its predecessor but core time regresses from 37.7s to 162.8s.
The truth table contains 409 zero-read barcodes and four inconsistent parent
totals; the report separates these from observed-barcode errors. Do not claim
universal accuracy or speed superiority from these experiments.

Read [the exact-search experiment](benchmark/four_condition_comparison/EXACT_SEARCH.md)
for the current implementation and validation. The remainder of this document
is the historical handoff for the parent branch; its speed figures and
statements that remaining accuracy cannot improve are not current conclusions.


Branch: **`feat/design-aware-scoring`** (9 commits ahead of `main`).
**Nothing is merged to `main` yet, deliberately.** The trade-offs in section A
should be settled first.

---

## A. Where each method actually wins and loses

Written from measurements in this document, not from the papers.

### barbac

**For.** The only method with no catastrophic regime — every competitor has one.
Ties the best available accuracy on the Johnson reference (FN 470, the lowest of
any method) and on anchored substitutions, and wins anchored+indel outright.
Fastest or near-fastest almost everywhere. Deterministic: identical input in any
row order gives identical output. Conserves every read. The only tool that also
does extraction, and it recovers 99.3% of reads where the previous in-house
pipeline recovered 0.64%. The only one offering both Hamming and Levenshtein.

**Against.** The default (`method="lv"`) is the wrong choice for the fixed-length
data that alignment-based extraction always produces — 117x slower for output
differing by 13 centroids in 174,000. Cost scales with *cluster count* rather
than input size, so diverse libraries are disproportionately expensive.
`merge_ratio = 20` cannot fire on unamplified libraries where nothing is 20x
anything, silently disabling error correction. The anchored-substitution parity
with Shepherd was obtained by implementing Shepherd's binomial criterion, so it
is "the same decision, 47x faster" and not an independent accuracy advantage.
Heaviest install of the four: an R package with a C++ core plus a conda
environment for the CLI stages.

### Shepherd

**For.** Genuinely best-in-class accuracy on substitution-only data; barbac had
to adopt its promotion rule to match it. Estimates its own error rate from the
data rather than taking a constant.

**Against.** Hamming only — it cannot represent an indel, and collapses when one
appears (118% and 257% wrong-split on our indel conditions). Its output depends
on input row order: permuting an identical file moves it between 471 and 478
false negatives on the Johnson data, so a single reported number is one draw.
Slow: 109s where barbac takes 23.6s, and 190s where barbac takes 4s. Loses reads
(88 of 25M on Johnson). Requires a fixed read length.

### Starcode

**For.** Robust — survives indels, which only it and barbac do. Deterministic.
A standalone C binary with no dependencies, which is the easiest thing on this
list to install. Marginally the best on random barcodes with indels.

**Against.** Four times worse on anchored designs (FN 319 against 79). On the
Johnson reference it is 1.6x worse on false negatives and 4x worse on false
positives than barbac. Scales badly: 3 minutes on 1.77M unique sequences,
over 17 minutes on 4.45M. Its default message-passing mode refuses to link any
pair unless the parent is 5x more abundant, which silently suppresses correct
merges on flat libraries -- we run `-s` (sphere) partly for that reason.

### Bartender

**For.** Frequently the fastest, and deterministic.

**Against.** Unusable with indels: 511%, 902% and 1486% wrong-split. Eight times
worse than barbac on false positives even on the Johnson data. Its interface
takes one row per *read* rather than per unique sequence, so a 25M-read
condition becomes a ~700MB file before it can start.

### The honest summary

barbac loses two categories by 6-12 barcodes and wins two by 210-240. Its
defensible claim is **"the only method that is never bad, and the fastest"** --
not "best in every category", which rests on margins smaller than the noise from
an arbitrary tie-breaking choice.

---

## B. CLI pipeline: verified end to end, two bugs fixed

Run on public data (`SRR9940679`, Johnson et al. 2019 lineage tracking,
BioProject `PRJNA559526`): FastQC, PEAR, minimap2, samtools, MultiQC,
`barbac_xtr` and `super_cluster2` all complete. 200,000 read pairs produced
61,006 merged reads, ~67% mapped, 17,884 full-length barcodes, clustered to
12,561 (Hamming, 0.57s) or 12,432 (LV, 5.04s).

Two bugs found and fixed (commit `6967399`):

1. **The pipeline could not run the tools barbac installs.**
   `configure_environment()` puts the tools in a conda environment;
   `run_cli_pipeline()` invoked them as bare commands through `PATH`. Every step
   failed with "command not found" -- and the failures were logged as warnings
   while the step still printed success, so the run continued two more stages
   and died with "No merged FASTQ files found", three steps from the cause.
2. **`check_barbac_tools()` had the mirror-image bug.** It searched only the
   conda environment (so system-wide tools were reported missing) and decided
   availability from whether `--version` succeeded -- which PEAR does not
   support, so an installed PEAR was reported absent.

**Undocumented trap worth adding to the docs:** a reference must not use `N` at
variable positions. A reference that was 44% `N` mapped **zero** reads; filling
those positions with a representative base mapped **67%** of the same reads at
the same coordinates. minimap2 cannot seed on `N`. The existing
`Reference_barcodes.fasta` works only because it is 24% `N`.

Also worth a look: PEAR merged only 61,006 of 200,000 pairs (30%). For a ~130bp
amplicon with 150bp mates that is low and may indicate a different insert size
than assumed.

---

## C. Immediate next step: a real time series

`PRJNA559526` contains **two 10-point time series** from the same evolution
experiment, which is the natural next test and the one that unlocks the only
remaining accuracy idea:

| series | library | timepoints | runs |
|---|---|---|---|
| C1 | `ILT_YPD_1` .. `ILT_YPD_10` | 10 | SRR9940651-60 |
| D1 | `ILT_YPA_1` .. `ILT_YPA_10` | 10 | SRR9940677-96 |

(`SRR9940679` already downloaded is D1 timepoint 4.)

**Plan.** Take one series -- D1 is smaller at ~53M reads total -- and subsample
each timepoint to ~500k read pairs, giving ~5M reads and roughly 400MB rather
than 4GB. Run the CLI pipeline per timepoint against the derived reference,
cluster each, then join on barcode to build trajectories.

**Why this is the right next step, in order of value:**

1. **It exercises `barbac_ts_area`**, the time-series half of the package, which
   nothing in this whole benchmarking effort has touched.
2. **It tests the one accuracy idea that adds information.** A barcode present
   at several timepoints is real; a sequencing error is not reproduced across
   independent library preparations. This is the only proposal that supplies
   evidence no competitor uses -- and unlike quality scores, it needs no change
   to extraction.
3. **It is a real-data result with a biological check.** Trajectories should be
   smooth and adaptive lineages should sweep. A clustering error shows up as a
   barcode that appears from nowhere or vanishes mid-experiment, which is a
   sanity test the single-sample benchmarks cannot provide.
4. **It is the discriminating real dataset we lack.** The ANC library is
   unamplified, so nothing is 20x anything and every method returns the input
   nearly unchanged. An evolved, amplified population has the abundance
   structure error correction actually depends on.

**Caveat to plan for:** the construct here was derived from the reads, not from
a published template. It should be re-derived per series and checked, since the
inline index length varies between libraries.

---
## 0. Latest result: barbac Hamming now exactly matches Shepherd

The user explicitly expects Hamming barbac to have at least Shepherd's accuracy
on Miloš's substitution-only algorithm. That expectation was reasonable, but
before this follow-up the implementations only shared Hamming distance; their
clustering decisions were not equivalent. This has now been corrected for the
two Shepherd decisions that mattered.

### What differed

Shepherd single-time-point clustering:

1. visits sequences by descending observed count;
2. chooses the nearest existing centroid (highest-count centroid breaks a
   distance tie);
3. merges every distance-1 neighbor;
4. merges singletons out to its learned `tau` (3 in this benchmark);
5. at distance 2/3, uses an exact binomial Bayes score with threshold `-4`;
6. preserves input order for equal counts.

barbac Hamming previously used a likelihood-scored parent, a count-ratio guard,
a conservative Poisson-style post-pass, and deterministic sequence tie order.
Calling both algorithms “Hamming” therefore did not make them identical.

### Exact diagnosis on anchored substitutions

Matched condition: 10,000 true barcodes, 1,000,000 reads, seed 42, lognormal
abundance, 0.5% substitutions/base, no indels, template
`NNNNNNNNATGCNNNNNNNNATCGTTAA`, truth-independent sequence tie order.

The old Hamming pass produced 9,975 roots. Its refinement correctly promoted
64 true barcodes and no false barcodes, reaching 10,039 centroids. However, 44
low-count errors had already founded roots because their correct parent did not
exist until that promotion pass. Every one of the 44 extra barbac centroids was
an error of a newly promoted true barcode:

- 37 were Hamming distance 1 (36 count-1, one count-2);
- six were distance 2 singletons;
- one was a distance 3 singleton.

Disabling refinement was tested and was decisively worse: FN increased from 83
to 147 while FP stayed 122. Lowering `merge_ratio` from 20 to 1 was also worse:
FN 84, FP 221, WS 212. Do not remove refinement or tune the global ratio.

### Implemented correction (`src/clustering.cpp`, build v11)

- In Hamming refinement only, replace the Poisson surrogate with Shepherd's
  exact single-time-point binomial Bayes decision (`bft = -4`). LV retains its
  previous rule and output.
- After promotion, revisit pre-existing roots only against newly promoted
  Hamming centroids. Absorb only Shepherd's unconditional safe cases: distance
  1, or a singleton within the configured distance. This avoids a naive second
  pass, which was shown to remove six real multi-read distance-3 barcodes.
- Add diagnostic counter `post_promotion_absorbed` and regression coverage.

Result:

| anchored substitutions | centroids | FN | FP | WS | algorithm time |
|---|---:|---:|---:|---:|---:|
| old barbac Hamming v10 | 10,039 | 83 | 122 | 113 | 4.03s |
| barbac Hamming v11 | **9,999** | **79** | **78** | **69** | **4.02s** |
| Shepherd | **9,999** | **79** | **78** | **69** | 190.37s |

The v11 barbac and Shepherd centroid sets are exactly equal, and every cluster
count is exactly equal. Thus barbac now reproduces Shepherd's complete result
about **47x faster** on this structured substitution benchmark. On random
substitutions, v11 remains FN 54 / FP 58 / WS 54 (no accuracy regression) and
runs in 0.65s versus Shepherd's 4.41s.

This does **not** prove universal equivalence. Shepherd estimates its error rate
(0.0049854 here), while barbac defaults to 0.005; Shepherd remains input-order
dependent on count ties, while barbac is deterministic. On the fair sequence-
ordered inputs used here those differences did not change the result. The dense
100k substitution benchmark has not yet been rerun with v11 because Hamming on
that 1.77-million-sequence input is slow; section 0.1 records the already-proven
tie-order diagnosis for its historical gap.

### Current four-condition snapshot

All conditions use 10,000 truths, 1,000,000 reads, seed 42, lognormal abundance,
0.5% substitutions/base, and sequence tie order. “Low indels” adds 0.5%
insertions and 0.5% deletions per base. Do not run high-indel or Nanopore
categories: the user explicitly excluded them.

| condition | barbac | Shepherd | Starcode | Bartender |
|---|---:|---:|---:|---:|
| random substitutions | Hamming **54/58/54; 0.65s**; LV 54/58/54; 1.19s | 54/58/54; 4.41s | **48/51/47**; 4.56s | 48/52/48; 1.63s |
| random substitutions + low indels | LV **170/350/171; 2.92s** | 103/4,924/4,890; 5.56s | **158/338/159; 30.83s** | 94/51,102/50,925; 4.33s |
| anchored substitutions | Hamming **79/78/69; 4.02s**; LV 80/78/69; 12.92s | **79/78/69; 190.37s** | 319/266/257; 6.80s | 105/108/100; 2.23s |
| anchored substitutions + low indels | LV **266/957/271; 38.99s** | 177/10,340/10,199; 163.73s | 476/1,180/494; 39.80s | 206/77,608/76,934; 10.65s |

Each cell is `FN / FP / WS; algorithm seconds`; lower is better. The code,
configuration, tool revisions, and compact CSV are under
`benchmark/four_condition_comparison/`. The orchestrator clean-builds barbac,
checks build v11, generates all data deterministically, and writes large raw
outputs under ignored `generated/`.

### Bottom line and next work for Claude

What succeeded:

- lossless LV A/C/G/T composition prefilter: 1.23x on `deep_sub_only`, 1.11x
  on `dense_sub_only`, byte-identical output;
- Hamming indel rescue on trace off-length reads: Johnson FP 647 -> 92;
- uninformative LV seed skipping: 4.8x on anchored designs, byte-identical;
- vectorized simulator: 4.8x faster with identical generated read multisets;
- Hamming Bayes promotion plus post-promotion cleanup: exact Shepherd output on
  anchored substitutions at ~47x its speed;
- fair sequence-based tie ordering exposed the old dense Shepherd advantage as
  simulator leakage rather than model accuracy.

What failed or should not be repeated:

- switching dense substitution data from LV to Hamming did not explain the old
  gap (FN 117 -> 113 only);
- removing refinement made anchored Hamming FN 147;
- `merge_ratio = 1` made anchored Hamming FP/WS explode to 221/212;
- tuning tie order or constants to the simulator's truth-first generation order
  is benchmark overfitting;
- marginal low-count truths are often absent or statistically indistinguishable
  from errors, so no clustering rule can recover all of them reliably.

Recommended next steps:

1. Claude owns the remaining speed work. Preserve byte-identical results while
   testing an absolute/sublinear LV posting-list cap, anchor-entropy-aware seeds,
   or a trie dynamic-programming search inspired by Starcode's `poucet` code.
2. Rerun dense substitution-only Hamming with v11 when runtime permits and check
   exact centroid sets against Shepherd under sequence tie order. Do not use the
   leaked generation-order score as a target.
3. For accuracy beyond Shepherd, use information Shepherd does not use: learned
   position/base/edit-specific error rates, actual alignment-path likelihoods,
   base qualities, joint replicate/time-point evidence, and reported uncertainty
   across deterministic tie seeds. Validate across multiple seeds and library
   densities; do not accept a change from one favorable benchmark.
4. User-facing method guidance should now be: random fixed-length substitutions
   -> Hamming; structured anchored substitutions -> Hamming v11 (now exact
   Shepherd accuracy and much faster); any real indels -> LV.

---

## 0.1 Earlier follow-up: dense gap resolved and lossless LV speedup

The apparent Shepherd accuracy advantage on `dense_sub_only` was a **simulator
tie-order leak**, not a Bayesian-model advantage. The simulator inserted each
true sequence into its `Counter` before its derived error variants and sorted
observations only by count. Shepherd preserves input order within count ties;
barbac deliberately uses a content-derived deterministic tie order.

Of the 40 true barcodes found by Shepherd but missed by barbac Hamming, 39 were
tied with the competing neighbor. The simulator put truth first in 39 of the 40
cases overall (38 of 39 tied cases), while barbac's lexicographic order put the
other sequence first in all 40. Randomizing the input before Shepherd's stable
abundance sort changed its result as follows:

| dense 100k result | FN | FP | WS |
|---|---:|---:|---:|
| Shepherd, simulator generation order | 73 | 134 | 73 |
| Shepherd, deterministic shuffled tie order | 112 | 173 | 112 |
| barbac Hamming, deterministic sequence order | 113 | 179 | 118 |

The dramatic 40-FN gap therefore collapses to one under a truth-independent tie
order. Do not tune barbac's merge guard or copy Shepherd's Bayes score to chase
the original number. `benchmark/indel_experiment/analyze_disagreements.py`
reproduces the set-level diagnosis. The simulator and runner now accept
`tie_order="sequence"` / `--tie-order sequence`; historical generation order
remains the simulator API default so old simulations still reproduce. The
comparison runner defaults to sequence order because its purpose is a fair
cross-method benchmark. It also defaults to `--barbac-method auto`, selecting
Hamming only when the simulated condition has zero insertions and deletions.

A lossless A/C/G/T-composition lower bound was also added before Levenshtein
verification. If the composition-vector L1 distance exceeds `2 * D`, edit
distance must exceed `D`. Measured results, with centroid CSVs byte-identical:

| dataset | old LV | composition prefilter | speedup | LV calls removed |
|---|---:|---:|---:|---:|
| `deep_sub_only` (176k sequences) | 3.84s | 3.12s | 1.23x | 52.1% |
| `dense_sub_only` (1.77m sequences) | 251.8s | 227.2s | 1.11x | 51.8% |

LV verbose output now recommends Hamming when every observed barcode has the
same length, explicitly conditional on indels and shifted alignments being
excluded. The default and clustering results are unchanged.

---

## 1. What is on the branch

| commit | what | verification |
|---|---|---|
| `d9fa30d` | Skip uninformative seeds in the LV fast path | centroid tables byte-identical on all 4 conditions; 4.8x faster on anchored designs |
| `9ab4220` | `tie_break` / `tie_seed` args + sensitivity analysis + figure | default output byte-identical |
| `cfb8763` | Rescue indel reads in Hamming mode | LV byte-identical; Hamming FP 647 -> 92 on Johnson 100k |
| `3d6f2a4` | `--abundance johnson` in `simulate.py` | default (`lognormal`) unchanged |
| `a45a935` | Vectorised the read simulator | identical read multiset at 0%, 0.5%, 2% indels; 4.8x faster |

**The invariant used throughout: any change to clustering must leave LV output
byte-identical unless it is deliberately changing behaviour.** Every commit above
was validated that way (`cmp` on the centroid CSVs). Please keep doing this — it
is the only reliable guard, and it is how a previous attempt at this work
(branch `fix/correct-parent-single-end`) was caught introducing an 11x slowdown
and a correctness bug.

---

## 2. Resolved historical problem — Shepherd on dense substitution data

The table below motivated the investigation. Section 0.1 now explains why its
largest gap is not an algorithmic accuracy difference.

| dataset | collisions* | indels | Shepherd FN | barbac FN | gap |
|---|---|---|---|---|---|
| `deep_sub_only` (10k barcodes) | 0.04% | no | 10 | 11 | 1 |
| Johnson et al. 100k (real published) | 0.52% | yes | 471 | 472 | 1 |
| `dense_sub_only` (100k barcodes) | 0.55% | **no** | **73** | **117** | **44** |

\* fraction of true barcodes lying within edit distance 3 of another true barcode.

The gap appears only when library density is high **and** there are no indels.
At low density the methods tie; on Johnson's data (high density, with indels)
they tie exactly — identical FP (86) and WS (61).

**Refuted hypothesis (do not repeat):** that LV's extra reach over-merges
near-neighbour true barcodes where Hamming would not. Tested directly —
`barbac method="hamming"` on `dense_sub_only` gives **FN 113**, essentially the
same as LV's 117, nowhere near Shepherd's 73. The distance metric is not the
cause.

**Refuted follow-up hypothesis:** the merge guard or Shepherd's Bayesian score
causes the 44-FN gap. Source inspection and set-level diagnostics show that the
gap consists almost entirely of equal-count orientation choices. Shepherd's
single-time-point implementation also merges distance-1 neighbors
unconditionally; its Bayes score only decides distance-2/3 cases below its
frequency threshold.

---

## 3. Confirmed findings (measurements, safe to build on)

**LV is the wrong default for fixed-length data — 86x cost for nothing.**
On real ANC amplicon data (178,213 barcodes extracted via `barbac_xtr`):

```
barbac Hamming   174,224 centroids   9.3s
barbac LV        174,211 centroids  798.3s     <- 13 centroids different, 86x slower
Shepherd         174,471 centroids   5.7s
```

`super_cluster2()` defaults to `method="lv"`. On alignment-extracted barcodes
(all one length, so no indels can be present) LV buys nothing. **A diagnostic
telling the user to switch is probably the highest-value small change left.**
Do not silently change the default — that would alter existing users' results.

**Residual clustering error is information-limited, not algorithmic.** Three
independent lines of evidence:
- 21 of 54 missed barcodes in `sub_only` never appear in the input at all.
- In the marginal set (count-1 absorbed at distance 1 by a small cluster) there
  are 32 true barcodes among 3,840 errors. A Poisson likelihood-ratio test
  separates them at best 32 TP / 845 FP. A "neighbourhood footprint" test
  (real barcodes should have error variants of their own) gives 5 TP / 159 FP —
  the footprint does not exist at these abundances (mean 0.2 neighbours).
- barbac and Shepherd produce *identical* FP and WS counts on Johnson's data;
  two unrelated algorithms converging on the same errors means the errors are
  in the data.
- The learned per-base error rate is **0.00489** against the hardcoded 0.005, so
  there is no miscalibration to exploit.

**Benchmark design matters more than the algorithms.** Which method wins is
determined by how densely packed the barcode library is:

| collisions | outcome |
|---|---|
| 0.04% | all four methods tie within 4 barcodes — the benchmark cannot discriminate |
| 0.52% | barbac and Shepherd lead; Starcode 63% worse on FN |
| 1.40% | barbac wins alone |

The originally published random conditions sat at 0.04%, which is why they
showed only noise. **Report collision density for every condition.**

**Tie-break choice is noise, and it is larger than the between-method
differences on shallow data.** Across 15 alternative tie orders barbac's FN
spans 45–53 on `sub_only` (median 48) while the shipped lexicographic order
gives 54; on the anchored design lexicographic is the *best* of 16 orders.
On Johnson's data the spread is 4 in 100,000. Do not tune the tie-break — it is
fitting to the benchmark. See `benchmark/tiebreak_sensitivity/`.

**Shepherd is order-dependent; barbac, Starcode and Bartender are not.**
Permuting identical input rows moves Shepherd's FN 44 -> 47. Worth stating in
the paper, neutrally.

---

## 4. Where barbac stands (all measured this session)

**Johnson et al. (2023) 100k benchmark** — the authoritative third-party test
(`~/Documents/Projects/Barcodes/barbac-benchmark/`, doi:10.5281/zenodo.7052124):

| method | Pearson R | FN | FP | WS | time |
|---|---|---|---|---|---|
| Shepherd | 1.000000 | 471 | 86 | 61 | 109.2s |
| barbac LV | 0.999998 | 472 | 86 | 61 | 53.1s |
| barbac Hamming | 1.000000 | 473 | 92 | 67 | **21.4s** |
| Bartender | 0.999995 | 494 | 725 | 700 | 28.4s |
| Starcode | 0.999992 | 771 | 353 | 328 | 163.7s |

barbac ties the best accuracy and is the fastest tool. The Hamming numbers are
*after* `cfb8763`; before it, Hamming gave FP 647 because indel-bearing reads
each founded their own cluster.

**Four simulated conditions at realistic depth** (`deep_*`, Johnson abundance,
250 reads/barcode): barbac wins the anchored+indel condition outright
(FN 9 vs Starcode 23, WS 9 vs 27), ties Shepherd on anchored+no-indel (both
FN 1, WS 0), is a close second on random+indel (Starcode WS 37 vs 71, barbac
2.8x faster), and sits in a four-way tie on random+no-indel.

Shepherd and Bartender collapse in **every** indel condition (118%, 257%, 902%,
1486% wrong-split). Only barbac and Starcode are usable there.

---

## 5. Speed: what is known and what is left

**Known bottleneck.** With `distance = 3` the pigeonhole principle forces the LV
seed index to seeds of `len/(D+1)` bases — ~6bp for 25bp barcodes. At 174,000
centroids a 6bp seed (4,096 values) holds ~42 centroids per bucket, and each
query probes 4 blocks x 7 shifts x several offsets. That is ~1,000 candidates
per read, and it is why LV cost scales with **cluster count**, not input size:
176k sequences producing 10k clusters takes 4.9s; 176k sequences producing
174k clusters takes 798s.

**Remaining ideas, in order of expected value:**
1. **Tighten the uninformative-seed threshold.** `query_specific_seed` currently
   skips buckets larger than `n_centroids / 8`. That *loosens* as the table grows
   — at 174k centroids it only skips buckets over 21,750. An absolute cap or a
   sublinear function is likely better. (Suspected but **not** measured.)
2. **Learn the design's constant positions.** Per-position base entropy would
   identify fixed anchors; seeding only on variable positions would make seeds
   far more discriminative on anchored designs. Bigger change, interacts with
   indels shifting positions.
3. **Trie-based LV search.** Starcode's exact “poucet” search performs dynamic
   programming over trie nodes and is the strongest source-backed candidate for
   replacing broad short-seed posting lists. This is a larger architectural
   experiment; preserve the byte-identical LV invariant while evaluating it.

---

## 6. Real data

`benchmark/102019_testBarcodes/` holds three PCR replicates of the ancestral
library (R1 only, 150bp, ~1.3M reads total). The same samples appear as
**already-mapped BAMs** in `~/Documents/Projects/Barcodes/sbw25-barcoding/`
(`data/bam/PCR{1,2,3}_ANC.assembled.fastq_sorted.bam`), which is the faster route.

**Extraction:** use `barbac_xtr(bam, ref_name="Reference_barcodes", start_pos=54,
end_pos=78)`, then drop barcodes containing `-` (partial coverage; 0.6% of
entries) and keep length 25. This recovers **99.3%** of reads. Do **not**
re-implement extraction by matching anchors in the FASTQ — the reads are
reverse-complemented, and requiring exact 25bp anchor spacing silently discards
every indel-bearing read (this mistake was made and cost several hours).

**What the real data can and cannot show.** All four methods return 40/46
Sanger-verified barcodes — which is every one present in the reads — and all
produce ~174,000 centroids from ~178,000 inputs. The library is *ancestral*,
so it is unamplified: every barcode sits at ~2.6 reads and nothing is 20x
anything, which means the default `merge_ratio = 20` cannot fire. **This dataset
cannot discriminate between clustering methods.** It is still valuable for the
extraction result and as a real-data sanity check.

For a real dataset that *can* discriminate, use Johnson et al. 2019 lineage
tracking, NCBI BioProject **PRJNA559526** — real reads from an amplified,
evolved population.

---

## 7. Reproducing the benchmarks

- Simulator: `benchmark/indel_experiment/simulate.py`, now with
  `--abundance {lognormal,johnson}`. Seeded; conditions regenerate deterministically.
- Four-way comparison harness: `benchmark/indel_experiment/run_experiment.py`
  (its `evaluate()` defines FN/FP/WS and is what every number here used).
- Tie-break sensitivity: `benchmark/tiebreak_sensitivity/` (R scripts + figure).
- External tools: `~/Documents/Projects/Barcodes/barbac-benchmark/tools/`.
- Bartender needs one row per *read*, not per unique sequence, so a 25M-read
  condition expands to a ~700MB file. Budget for it or skip it — its result on
  indel data is not in doubt.

---

## 8. Known liabilities in the repo (not yet fixed)

1. **`README.md` benchmark tables do not reproduce.** The 10k tables were
   generated before commit `af933ec` (the tie-break fix) and current code gives
   different numbers (e.g. `sub_only` FN 44 published vs 54 actual). The 100k
   parity table *does* still reproduce.
2. **The Johnson dataset is cited without its DOI** and described as "the
   Johnson et al. (2023) reference dataset", which reads as though it were real
   sequencing reads. It is their *simulated* benchmark; cite
   doi:10.5281/zenodo.7052124 and say so.
3. Shepherd is attributed loosely — it is Tavakolian & Frazão et al.,
   *Bioinformatics* 2022, not part of Johnson et al.
