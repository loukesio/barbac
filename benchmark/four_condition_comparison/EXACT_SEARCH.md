# Exact search and support-based barcode ordering

This experiment starts from `0cbe94c` on `feat/design-aware-scoring` and is
implemented on `feat/exact-search-clustering`. It keeps the four existing
conditions, distance three, merge ratio 20, error rate 0.005 and the original
simulator. No truth labels enter clustering or sequence ordering.

The historical handoff correctly identified input-order leakage and the
limitations of fixed-length Hamming distance. It also left two useful avenues:
better search partitions and information beyond an arbitrary sequence tie.
This branch implements both, without altering the merge thresholds.

## What changed

**Information-balanced partitions.** Contiguous equal-sized blocks can fall
entirely inside a constant anchor. They then match nearly every centroid.
The new index estimates per-position collision information from observed base
counts. For Hamming it spreads informative positions across disjoint blocks.
For LV a small dynamic program chooses contiguous boundaries that minimise the
sum of estimated block collision probabilities. Correlations between positions
may weaken the cost estimate; they cannot invalidate candidate recall.

The exact-search guarantee is separate from that estimate: with D edits and
D+1 disjoint query blocks, at least one block survives unchanged. LV indexes
all permitted shifts from -D to D for each observed query length. Buckets have
separate block/length namespaces and full 64-bit sequence keys. Unsupported
LV strings fall back to scanning, including ambiguous characters and long
sequences. The accelerated Hamming domain remains A/C/G/T, at most 32 bases.

**A provable fast path.** LV first searches the one-edit index. It skips the
wider index only when the selected parent's score exceeds an upper bound for
every unseen parent: distance two, the maximum input abundance and the minimum
comparison length. Otherwise it searches the entire configured radius. The old
code stopped at any absorbable Hamming or long-seed match, so it could miss a
better indel-derived parent. Equal scores now resolve by centroid creation order
in both the indexed and exhaustive LV backends.

**Correct distance evaluation.** Hamming distance is an upper bound on LV,
not generally its value. Three mismatches can cost two edits. Packed Hamming
is now accepted as the exact LV value only at distance at most two. The long
sequence banded implementation also now initialises the right-hand band
boundary, which previously could read uninitialised stack values.

**Refinement cost.** Member reassignment visits only actually promoted members.
Previously it iterated over the entire cluster for each member, even when there
were no promotions. This removes quadratic work without changing the rule.

**Evidence for count ties.** `tie_break = "support"` ranks equal-count sequences
by the summed counts of their one-edit neighbours whose counts are no larger
than the sequence itself. More abundant neighbours are excluded so that a
nearby large lineage does not by itself establish a low-count sequence as real.
Remaining ties use sequence order. Hamming support uses Hamming distance; LV
support uses edit distance. This is optional; `"sequence"` remains the default,
and the existing `"hash"` sensitivity option remains available.

This draws on the useful idea of neighbourhood evidence seen in Starcode's
sphere-order comparator, but does not copy its complete distance-three graph
or unrestricted neighbourhood count. See the original
[Starcode paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4765884/) and
[source](https://github.com/gui11aume/starcode/blob/8987b2eb7558bdb6e8a44823881e4e2f76f47535/src/starcode.c).
The Hamming promotion rule inherited from the baseline still uses the
[Shepherd model](https://pmc.ncbi.nlm.nih.gov/articles/PMC9344852/).

## Using the change

```r
# Substitution-only, fixed-length data:
result <- super_cluster2(counts, method = "hamming", tie_break = "support")

# Data with insertions, deletions, or shifted alignments:
result <- super_cluster2(counts, method = "lv", tie_break = "support")
```

The index learns search partitions whether or not `use_design` is enabled.
`kmer_size` remains accepted for compatibility but no longer sets the partitions.
The public result columns are unchanged. The experimental design-based merge
rule and the default error/abundance parameters remain unchanged.

## Validation design

The development seed is 42. Seeds 43 and 44 are independent validation draws,
with 10,000 true barcodes and one million reads in every category. All tools see
identical sequence-ordered inputs. Comparators are Shepherd revision
`83ea670b929d4c3370e9f162bfdc34e6893cf5cc` and Starcode revision
`8987b2eb7558bdb6e8a44823881e4e2f76f47535`, sphere clustering at distance three.
Both barbac distances are measured in all categories, including Hamming's
unsuitable indel conditions.

`compare_revisions.py` checkpoints results after each run and records simulation
configuration, input checksums, source checksums, installed-package checksums,
competitor revisions, native build IDs, and environment information. Repeated
outputs must be identical and every barbac run must conserve all input reads.
`time_revision.R` provides separate repeated measurements in a loaded R session;
these exclude startup and must not be described as end-to-end runtime.

FN counts true barcode strings missing from the reported centroids. FP counts
reported centroid strings absent from truth. WS is the subset of FP within
three edits of a true barcode. Centroid F1 is `2 TP / (2 TP + FN + FP)` and
balances missed barcodes against spurious ones. These are centroid metrics,
not read-assignment accuracy. Truth strings never observed exactly are reported
separately: a method restricted to observed centroids cannot recover them.

The regression suite checks full-scan equivalence over multiple lengths and
radii, explicit better-parent and cyclic-shift counterexamples, ambiguous and
long-string fallback, read conservation, support against an independent distance
matrix, and invariance to input order. Ordinary package tests are also run.

## Measured accuracy: three seeds, four categories

The support rule has the best **mean centroid F1** in all four categories across
seeds 42, 43 and 44. It improves FN and FP relative to the old barbac in every
seed/category under the indicated distance method. This is an observed result
on these simulations, not a confidence bound or a universal ranking.

| Condition | barbac distance | Old barbac F1 | New support F1 | Shepherd F1 | Starcode F1 |
|---|---|---:|---:|---:|---:|
| Random substitutions | hamming | 99.502% | **99.538%** | 99.502% | 99.475% |
| Random + indels | lv | 97.405% | **97.730%** | 79.648% | 97.456% |
| Anchored substitutions | hamming | 99.162% | **99.249%** | 99.162% | 97.039% |
| Anchored + indels | lv | 94.120% | **94.799%** | 65.188% | 91.962% |

LV also improves on both competitors' mean F1 on the substitution-only data:
99.538% for random barcodes and 99.244% for anchored barcodes. Hamming remains
the faster choice for known substitution-only data.

The historical seed-42 gaps shrink or reverse: random substitutions improve
from FN/FP 54/58 to 49/53, random indels from 170/350 to 134/314, anchored
substitutions (Hamming) from 79/78 to 68/67, and anchored indels from 266/957 to
200/886. **Starcode still narrowly wins seed-42 random substitutions**, 48/51;
the three-seed mean is better for barbac because seeds 43 and 44 also matter.

The exact search with the original sequence tie rule retains the baseline FN
in all 24 seed/category/distance comparisons. Its FP is also unchanged except
one fewer false centroid on seed-43 anchored indels. Cluster counts/assignments
are allowed to change where exact parent scoring corrects the old heuristic;
this branch does not claim universal byte-identical output.

All 72 barbac runs conserve all one million input reads. Shepherd loses roughly
24,000 reads per random-indel input and 45,000 per anchored-indel input in this
runner. Its lower FN on indels must therefore be read alongside its thousands
of FP and its read loss. Pure Hamming barbac retains roughly 55,000 random-indel
FP and 82,000 anchored-indel FP even with support: these are failures, and are
included in the full results rather than omitted.

The checked [96 per-run rows](summary_2026-09-08.csv),
[aggregates](aggregate_2026-09-08.csv), and
[provenance manifest](versions_2026-09-08.json) include both distances in all
categories. The barbac and competitor campaigns ran separately and concurrently;
all 12 pairs of independently generated input files were verified byte-identical.
The manifest identifies the timing limits and stores input/output checksums.

The largest remaining anchored-indel FP component is outside the current
radius: mean FP is 882.3, while mean WS is 194.3, leaving **688 false centroids
farther than three edits from every true barcode**. That motivates evaluating a
carefully controlled tail-error rescue, rather than changing a tie constant.

## Measured speed and its limits

Median of three sequential repetitions on seed 42, in a loaded R session.
Each timed call includes CSV input and the R wrapper. The accuracy and
competitor campaigns had finished before these runs; part of the baseline run
overlapped package checks. These are engineering measurements, not controlled
hardware-isolation experiments.

| Condition | Distance | Old sequence rule | New sequence rule | New support rule |
|---|---|---:|---:|---:|
| Random substitutions | hamming | 0.245s | 0.236s | 0.370s |
| Random substitutions | lv | 0.920s | 0.497s | 0.846s |
| Random + indels | lv | 3.124s | 1.920s | 3.887s |
| Anchored substitutions | hamming | 3.581s | 0.460s | 0.656s |
| Anchored substitutions | lv | 12.702s | 0.995s | 1.683s |
| Anchored + indels | lv | 38.593s | 4.904s | 8.619s |

With the accuracy-improving support option, anchored Hamming is **5.5x faster**,
anchored substitution LV **7.5x faster**, and anchored-indel LV **4.5x faster**
than the baseline. The exact index alone gives larger speedups. Support has a
real cost: random-substitution Hamming rises from 0.245s to 0.370s, and random
indel LV from 3.124s to 3.887s. It exchanges some of the indexing speed gain for
better centroid recovery. All [72 timing observations](timing_2026-09-08.csv)
are retained, including unsuitable Hamming-indel runs.

**End-to-end speed is not best in every category.** In the separate three-seed
campaign, support Hamming takes a median 5.76s as a fresh process on random
substitutions, versus Starcode 4.64s and Shepherd 5.21s. Support LV on random
indels takes 9.42s including startup, versus Shepherd 5.43s (with much worse
accuracy) and Starcode 28.63s. R startup and package loading account for much of
the difference. The full CSV retains both timings; warm clustering time must
not be compared to another tool's end-to-end time as though they were identical.

## Validation assessment: share with caveats

`R CMD build --no-manual` and `R CMD check --no-manual` completed with tests,
examples, documentation checks and rebuilt vignettes passing. The
[package check log](check_2026-09-08.log) records one warning from the installed
R header's unsupported `-Wfixed-enum-extension` warning group under the local
Apple compiler. Repository index connectivity was unavailable in the sandbox;
installed dependency checks completed. There were no package-check errors.

The headline set-based metrics, F1 arithmetic, input identity and read totals
were independently recomputed for all 96 rows. The 72 repeated timing records
were checked for completeness, conserved reads and identical repeated outputs.
Remaining caveats are the three-seed simulation scope, centroid rather than
read-assignment accuracy, unsuitable Hamming-indel behavior, and the distinction
between loaded-session and end-to-end timing. No best-on-every-seed or
best-end-to-end-in-every-category claim is supported.

## Further experiments with the highest value

1. Learn separate substitution, insertion and deletion likelihoods from
   high-confidence error clouds. A scalar edit penalty treats distinct error
   paths equally and is especially crude in repeats. Fit on one subset of
   abundant lineages and test on withheld lineages before changing merges.
2. Record per-base quality during extraction and combine it with alignment
   likelihood. Equal-count, weakly supported barcodes often cannot be separated
   from strings alone; quality adds independent evidence without changing the
   synthetic tie order. This requires quality-bearing inputs and a benchmark
   that actually simulates or retains quality.
3. Assess a conservative, separately reported rescue pass for errors beyond the
   chosen radius, using a supplied barcode design when available. Increasing
   the global distance also increases true-lineage collisions, so it needs
   explicit false-merge measurements and validation on other library densities.

None of these unimplemented experiments is counted as an achieved improvement.
A uniform-rate simulator, one abundance distribution and one anchored design
cannot establish universal superiority. Pure Hamming cannot become a general
indel metric; use LV when the data requires it.
