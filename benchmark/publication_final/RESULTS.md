# Final publication benchmark

[One-page table](benchmark_table.pdf) · [Statistical supplement](benchmark_supplement.pdf) · [Full per-input results](publication_results.csv) · [Protocol](README.md)

60 independent libraries per simulated design; all 726 registered tool cells attempted, 10 failed. 6 of eight specified LV-versus-competitor F1 contrasts meet both lower-bound criteria; 2 are unavailable. Only full 60-pair contrasts are tested under the disclosed [failure-reporting addendum](FAILURE_REPORTING.md). Milo is a fixed reference, not an independent replication test.

## Random N20: substitutions + repeat-dependent indels

| Method | FN | FP | F1 % (SD) | Workflow seconds |
|---|---:|---:|---:|---:|
| barbac Hamming | 124.92 | 18.30 | 99.28004 (0.0758) | 1.40 |
| barbac LV + Poisson | 124.90 | 18.25 | 99.28038 (0.0752) | 1.76 |
| Shepherd [55/60] | n/a | n/a | n/a | n/a |
| Starcode sphere | 129.72 | 22.32 | 99.23569 (0.0797) | 3.29 |
| Starcode MP | 118.18 | 142.98 | 98.69582 (0.1112) | 3.22 |
| Bartender | 125.17 | 48.03 | 99.13063 (0.0921) | 1.89 |

## Anchored N20 + AA/TT: substitutions + repeat-dependent indels

| Method | FN | FP | F1 % (SD) | Workflow seconds |
|---|---:|---:|---:|---:|
| barbac Hamming | 127.70 | 24.90 | 99.23302 (0.0736) | 1.52 |
| barbac LV + Poisson | 127.67 | 24.08 | 99.23726 (0.0741) | 2.05 |
| Shepherd [55/60] | n/a | n/a | n/a | n/a |
| Starcode sphere | 132.63 | 27.65 | 99.19431 (0.0829) | 3.77 |
| Starcode MP | 119.95 | 191.83 | 98.44674 (0.1384) | 3.67 |
| Bartender | 127.73 | 137.97 | 98.67218 (0.0786) | 2.31 |

## Milo / Johnson published reference simulation

| Method | FN | FP | F1 % (SD) | Workflow seconds |
|---|---:|---:|---:|---:|
| barbac Hamming | 469 | 87 | 99.72147 | 17.50 |
| barbac LV + Poisson | 471 | 83 | 99.72246 | 33.65 |
| Shepherd | 470 | 85 | 99.72196 | 120.12 |
| Starcode sphere | 771 | 353 | 99.43682 | 170.68 |
| Starcode MP | 532 | 560 | 99.45408 | 193.39 |
| Bartender | 496 | 795 | 99.35546 | 41.77 |

## Accuracy contrasts

| Design | Competitor | Mean F1 difference (points) | Simultaneous lower bound | Bootstrap lower bound | Supported |
|---|---|---:|---:|---:|---|
| random_mixed | Shepherd | n/a | n/a | n/a | n/a |
| random_mixed | Starcode sphere | +0.04468 | +0.03322 | +0.03372 | Yes |
| random_mixed | Starcode MP | +0.58456 | +0.55942 | +0.56049 | Yes |
| random_mixed | Bartender | +0.14975 | +0.13852 | +0.13915 | Yes |
| anchored_mixed | Shepherd | n/a | n/a | n/a | n/a |
| anchored_mixed | Starcode sphere | +0.04295 | +0.03084 | +0.03138 | Yes |
| anchored_mixed | Starcode MP | +0.79053 | +0.75721 | +0.75913 | Yes |
| anchored_mixed | Bartender | +0.56509 | +0.54439 | +0.54550 | Yes |

## Interpretation and scope

Bold/shading identifies numerical leaders among complete rows, not significance. 10 failed tool calls are retained. An incomplete row shows successes/planned in brackets; n/a means the full-design result is unavailable. Conditional summaries and all eight accuracy contrasts appear in the supplement.
FN: true identities absent from inferred centroids. FP: inferred identities absent from truth. Exact-identity F1 includes zero-read true barcodes. Read assignment, positive-read truth and abundance metrics are reported separately.
New simulations: 60 independent seeds, 10,000 true identities and one million expected reads per library; paired designs share variable identities and parent counts. Substitutions: 0.4% per base. Indels: unscaled archived homopolymer-specific rates.
The calibrated table ends at repeat length 13. Further indel recursion is stopped for affected reads beyond that support; reads, substitutions and truth are retained. Boundary use: 249 reads across 16 of 120 libraries. This is a declared finite-support model, not an empirical estimate beyond the table.
Milo: unchanged deposited simulation, 100,000 true barcodes and 24,996,128 reads; used during development. Its single-reference scores do not establish independent-library statistical superiority. Its timing is one observation in this final campaign.
All methods use distance 3 and one thread. Barbac: support ordering, ratio 20, configured error proxy 0.005, design scoring off, LV Poisson on. Hamming retains rare-indel rescue. Shepherd: nominal length, Bayes threshold 4. Starcode MP ratio 5. Bartender: seed length 5, step 1, z=5, cutoff 1.
Timing: serial fresh-worker elapsed time including startup, required conversion and centroid/member exports; staging, scoring and hashing excluded. The 28 calls overlapping a system sleep/partial-wake interval were repeated once with sleep prevented; mappings had to match exactly. Original and repaired records are retained. Other cells ran once.
The publication scope was selected after development inspection. Code, settings, seeds and test methods were frozen before final generation. A disclosed reporting addendum retains tool failures and tests only fully observed 60-pair contrasts, keeping the eight-comparison correction. All registered libraries remain. Development controls are preserved.

The [timing repair](TIMING_REPAIR.md) retains both original and selected measurements. Original accuracy and timing receipts remain in [all_results.csv](all_results.csv).

## Matched successful-only context

These are descriptive comparisons on the same inputs where the competitor completed. They do not replace full-design tests, and failed inputs remain retained.

| Design | Completed inputs | LV F1 % | Competitor | Competitor F1 % |
|---|---:|---:|---|---:|
| random_mixed | 55 | 99.29302 | Shepherd | 99.05236 |
| anchored_mixed | 55 | 99.25247 | Shepherd | 99.23909 |
