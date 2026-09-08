**3.3 Comparison with established error-correction methods**

We compared barbac Hamming and Levenshtein (LV) clustering with Shepherd, Starcode
sphere, Starcode default message passing, and Bartender using four simulation conditions
and the Johnson et al. (2023) reference dataset, referred to here as Milos (Table 2).
The four conditions crossed fully random 20-base barcodes with anchored 28-base barcodes
containing 16 variable positions, and substitutions alone with substitutions plus
indels. Each condition used 10,000 true barcodes and one million reads with lognormal
abundances (sigma = 1.5), independently generated with seeds 42, 43 and 44. The
substitution probability was 0.005 per base; indel conditions additionally used
insertion and deletion probabilities of 0.005 per opportunity each. These simulations
evaluate specified error regimes rather than representing every sequencing platform. The
Milos input contains 1,544,850 unique observed sequences and 24,996,128 reads for
100,000 listed true barcode sequences.

All comparisons used maximum distance three. Both barbac modes used `super_cluster2()`,
native build v13, merge ratio 20, configured error rate 0.005, support ordering for
equal-count sequences, and the design option disabled. LV additionally enabled
`indel_model = "poisson"`, an experimental exception for repeated-base single indels
whose abundance is consistent with an expected-error model. This option is disabled by
default; the table explicitly evaluates the enabled configuration. The error rate was
supplied, not fitted to true labels. Shepherd and Bartender used their remaining default
parameters, with Shepherd estimating its substitution error rate automatically. Starcode
sphere and default message passing were tested separately. Timed Starcode and Bartender
runs each requested one thread. Input sequences were ordered by decreasing observed
count and then sequence, without using parent labels. The v13 LV search optimization
excludes candidate parents whose best possible likelihood score cannot improve the
current assignment; the indexed results were checked against existing outputs and
full-scan tests. Hamming refinement incorporates the binomial criterion described for
Shepherd (Tavakolian et al., 2022).

For this comparison, a true positive (TP) is an output centroid string present in the
true barcode set; FN counts missing true strings and FP counts extra centroid strings.
Centroid F1 is `2 TP / (2 TP + FN + FP)`. All listed truth sequences remain in the
denominator, including those whose exact sequence does not occur in the noisy input. For
the Milos simulation this includes 409 truth barcodes with zero reads. Input identity,
output uniqueness and count reconciliation were checked before scoring. All barbac
outputs conserved the supplied reads. Per-read parent labels were not retained for the
four smaller simulations, so no read-assignment accuracy is inferred from their centroid
F1 values.

**Table 2. Latest barcode clustering performance at distance three.** R-S: random
barcodes with substitutions; R-I: random barcodes with substitutions and indels; A-S and
A-I: the corresponding anchored designs. Milos denotes the Johnson et al. (2023)
reference simulation. FN and FP are absolute barcode counts; for R-S, R-I, A-S and A-I,
counts and centroid F1 are means across seeds 42, 43 and 44 (10,000 true barcodes and
one million reads per seed). The Milos row for each method uses one fixed dataset of
100,000 true barcodes and 24,996,128 reads. Bold F1 values identify the highest score
within a dataset. F1 measures exact centroid recovery and differs from the cluster-label
F1 used in Section 3.1. Time is one serial workflow observation in seconds on the same
development machine: seed 42 for the four smaller simulations and the full dataset for
Milos. It includes program startup, required format conversion, clustering, and
centroid/member exports; shared input staging and scoring are excluded. These timings
have no confidence intervals. Both barbac modes use native v13 and support ordering; LV
additionally enables the experimental Poisson indel option. Starcode sphere and default
message passing (MP) are reported separately. All methods use distance three, with its
method-specific Hamming or Levenshtein interpretation.

| Dataset | Method | FN | FP | F1 (%) | Time (s) |
|---|---|---:|---:|---:|---:|
| R-S | barbac Hamming | 44.7 | 47.7 | 99.538 | 4.98 |
| R-S | barbac LV + Poisson | 44.7 | 47.7 | 99.538 | 4.14 |
| R-S | Shepherd | 48.3 | 51.3 | 99.502 | 3.47 |
| R-S | Starcode sphere | 51.3 | 53.7 | 99.475 | 2.94 |
| R-S | Starcode MP | 25.7 | 430.0 | 97.767 | 3.12 |
| R-S | Bartender | 44.3 | 46.3 | **99.547** | 1.57 |
| R-I | barbac Hamming | 98.7 | 55,318.7 | 26.326 | 4.78 |
| R-I | barbac LV + Poisson | 134.7 | 322.3 | **97.736** | 5.96 |
| R-I | Shepherd | 104.0 | 4,953.3 | 79.648 | 3.74 |
| R-I | Starcode sphere | 163.3 | 350.3 | 97.456 | 17.52 |
| R-I | Starcode MP | 72.7 | 1,880.7 | 91.044 | 18.84 |
| R-I | Bartender | 98.3 | 51,037.7 | 27.916 | 3.05 |
| A-S | barbac Hamming | 72.3 | 78.0 | **99.249** | 4.05 |
| A-S | barbac LV + Poisson | 73.0 | 78.3 | 99.244 | 4.63 |
| A-S | Shepherd | 81.0 | 86.7 | 99.162 | 117.94 |
| A-S | Starcode sphere | 320.0 | 270.7 | 97.039 | 4.44 |
| A-S | Starcode MP | 146.0 | 621.7 | 96.251 | 5.32 |
| A-S | Bartender | 115.0 | 118.7 | 98.832 | 2.17 |
| A-I | barbac Hamming | 159.7 | 81,554.3 | 19.410 | 5.77 |
| A-I | barbac LV + Poisson | 193.7 | 882.0 | **94.801** | 8.50 |
| A-I | Shepherd | 170.7 | 10,327.7 | 65.188 | 105.08 |
| A-I | Starcode sphere | 476.0 | 1,189.0 | 91.962 | 26.24 |
| A-I | Starcode MP | 223.0 | 3,272.0 | 84.837 | 26.51 |
| A-I | Bartender | 207.3 | 77,390.7 | 20.153 | 7.90 |
| Milos | barbac Hamming | 469 | 87 | 99.72147 | 20.81 |
| Milos | barbac LV + Poisson | 471 | 83 | **99.72246** | 43.40 |
| Milos | Shepherd | 470 | 85 | 99.72196 | 110.95 |
| Milos | Starcode sphere | 771 | 353 | 99.43682 | 155.30 |
| Milos | Starcode MP | 532 | 560 | 99.45408 | 150.22 |
| Milos | Bartender | 494 | 725 | 99.39120 | 36.52 |

The highest-scoring barbac configuration achieved the highest mean centroid F1 in 3 of
the four simulation conditions, and LV achieved the highest centroid F1 on Milos (Table
2). Bartender narrowly led the random substitution-only condition. On random
substitution-only data, Hamming and LV had identical centroid recovery, while Hamming
slightly outperformed LV on anchored substitution-only data. The main distinction
emerged in indel-containing data: LV produced far fewer false clusters than Hamming,
Shepherd, and Bartender. For random indel data, LV yielded mean FN 134.7 and FP 322.3,
with F1 97.736%. For anchored indel data, the corresponding values were FN 193.7, FP
882.0, and F1 94.801%. Hundreds of false clusters remained in the latter condition, so
relative superiority does not imply error-free reconstruction. The Poisson exception
removed five false clusters across these 12 smaller datasets without changing FN; most
of the observed difference from competitors reflects the broader LV clustering strategy
and support ordering.

On Milos, LV returned FN 471 and FP 83, giving centroid F1 99.72246%, compared with
Shepherd's FN 470, FP 85, and F1 99.72196%. This F1 difference is very small and is
reported descriptively. With the supplied read-parent labels, LV made 131 wrong
assignments and left no reads unassigned; Shepherd made 129 wrong assignments and left
88 reads unassigned. Thus the combined wrong-or-unassigned count was 131 for LV and 217
for Shepherd. The reference labels contain four discrepancies in per-parent totals,
spanning six reads in absolute difference; these source inconsistencies were retained
and audited.

**3.4 Runtime and operational cost**

Runtime depended on dataset size, sequence design and the reported boundary. On Milos,
the LV workflow took 43.4 s, compared with 111.0 s for Shepherd (2.56-fold shorter).
Hamming took 20.8 s and was the fastest tested workflow on that dataset. On the smaller
simulations, interpreter/package startup contributed several seconds to barbac
workflows, and native competitors were sometimes faster end to end despite poorer
recovery. The table therefore reports complete measured workflows rather than comparing
barbac core time against another method's total runtime. Times for the four simulations
are fresh serial seed-42 observations; accuracy additionally includes seeds 43 and 44.
Previously recorded peer accuracy outputs were reused only after verifying identical
inputs and output hashes. Earlier overlapping benchmark timings were excluded. Timing
repetitions sufficient for confidence intervals were not collected.

**3.5 Scope and limitations of the comparison**

The results support an accuracy advantage for the tested LV configuration in the
indel-containing simulations, with Hamming useful for substitution-only inputs. They do
not establish universal superiority across libraries or sequencing platforms. The Milos
dataset has only 1,001 reads whose lengths differ from their supplied parents, whereas
the four-condition experiment explicitly includes broader mixed indel errors. Shepherd
includes separate correction of simple single insertions and deletions around a supplied
barcode length; its poorer performance in mixed-error simulations should not be
described as complete absence of indel support. The Poisson option can merge real length
variants with error-like counts, and its configured error rate is not a learned
platform-specific indel rate. Independent labeled experimental controls, additional
abundance distributions and repeated timings are needed before broader accuracy or speed
claims.

