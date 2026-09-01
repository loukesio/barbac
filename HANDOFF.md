# Handoff: barbac clustering performance and benchmarking

Branch: **`perf/clustering-improvements`** (5 commits ahead of `main`, working tree clean)

Everything below was measured on 2026-09-01 on the author's machine. Where a
claim is a measurement it says so; where it is a hypothesis it says that too.
Two hypotheses in here were tested and **refuted** — please read those before
re-deriving them.

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

## 2. The open problem — this is the interesting one

**barbac loses to Shepherd on dense, substitution-only libraries, and nobody
knows why.**

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

**Where to look instead:** the merge guard (`effective_merge_ratio` /
`effective_count_floor` in `src/clustering.cpp`) or the likelihood best-parent
scoring, specifically their behaviour when two *true* barcodes are near
neighbours. That is the only condition that distinguishes the three rows above.
Shepherd uses a Bayesian abundance test rather than a count-ratio rule; the
difference in how the two decide "is this a new barcode or an error of that one"
is the likeliest source.

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

**Untried ideas, in order of expected value:**
1. **Diagnostic when all input lengths are equal** — recommend `method="hamming"`.
   The length histogram is already computed in `clustering.cpp` (added in
   `cfb8763` for the Hamming indel warning). Cheap, safe, ~80x for affected users.
2. **Composition prefilter** — compare A/C/G/T counts before the Myers DP; if the
   count difference exceeds 2*D the edit distance must too. Lossless, a few
   instructions instead of a full DP. Should cut roughly half the distance
   computations.
3. **Tighten the uninformative-seed threshold.** `query_specific_seed` currently
   skips buckets larger than `n_centroids / 8`. That *loosens* as the table grows
   — at 174k centroids it only skips buckets over 21,750. An absolute cap or a
   sublinear function is likely better. (Suspected but **not** measured.)
4. **Learn the design's constant positions.** Per-position base entropy would
   identify fixed anchors; seeding only on variable positions would make seeds
   far more discriminative on anchored designs. Bigger change, interacts with
   indels shifting positions.

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
