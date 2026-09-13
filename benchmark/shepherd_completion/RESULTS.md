# Complete comparison with Shepherd sensitivity configuration

The original confirmatory benchmark is unchanged. Shepherd uses the separately registered post hoc configuration explained in [README.md](README.md). All other scores and accepted times are reused.

[PDF table](benchmark_table.pdf) · [CSV table](benchmark_table.csv) · [execution receipts](results.json)

## Random N20: substitutions + repeat-dependent indels

| Method | FN | FP | F1 % | Workflow s | Correct reads % | Wrong reads | Unassigned |
|---|---:|---:|---:|---:|---:|---:|---:|
| barbac Hamming | 124.92 | 18.30 | 99.28004 | 1.40 | 99.997077 | 29.23 | 0.00 |
| barbac LV + Poisson | 124.90 | 18.25 | 99.28038 | 1.76 | 99.997083 | 29.17 | 0.00 |
| Shepherd* | 126.53 | 19.83 | 99.26421 | 4.67 | 99.996495 | 32.00 | 3.05 |
| Starcode sphere | 129.72 | 22.32 | 99.23569 | 3.29 | 99.982520 | 174.78 | 0.00 |
| Starcode MP | 118.18 | 142.98 | 98.69582 | 3.22 | 99.983934 | 160.65 | 0.00 |
| Bartender | 125.17 | 48.03 | 99.13063 | 1.89 | 99.930972 | 690.32 | 0.00 |

Simulations: means over all 60 libraries; time is median. Milo: one observation. Shepherd was timed in a later session.

## Anchored N20 + AA/TT: substitutions + repeat-dependent indels

| Method | FN | FP | F1 % | Workflow s | Correct reads % | Wrong reads | Unassigned |
|---|---:|---:|---:|---:|---:|---:|---:|
| barbac Hamming | 127.70 | 24.90 | 99.23302 | 1.52 | 99.993137 | 68.62 | 0.00 |
| barbac LV + Poisson | 127.67 | 24.08 | 99.23726 | 2.05 | 99.994257 | 57.42 | 0.00 |
| Shepherd* | 129.65 | 25.05 | 99.22239 | 13.45 | 99.993818 | 39.37 | 22.45 |
| Starcode sphere | 132.63 | 27.65 | 99.19431 | 3.77 | 99.986583 | 134.17 | 0.00 |
| Starcode MP | 119.95 | 191.83 | 98.44674 | 3.67 | 99.979174 | 208.25 | 0.00 |
| Bartender | 127.73 | 137.97 | 98.67218 | 2.31 | 99.906649 | 933.50 | 0.00 |

Simulations: means over all 60 libraries; time is median. Milo: one observation. Shepherd was timed in a later session.

## Milo / Johnson published reference simulation

| Method | FN | FP | F1 % | Workflow s | Correct reads % | Wrong reads | Unassigned |
|---|---:|---:|---:|---:|---:|---:|---:|
| barbac Hamming | 469.00 | 87.00 | 99.72147 | 17.50 | 99.999088 | 228.00 | 0.00 |
| barbac LV + Poisson | 471.00 | 83.00 | 99.72246 | 33.65 | 99.999476 | 131.00 | 0.00 |
| Shepherd* | 470.00 | 85.00 | 99.72196 | 165.23 | 99.999132 | 129.00 | 88.00 |
| Starcode sphere | 771.00 | 353.00 | 99.43682 | 170.68 | 99.847880 | 38024.00 | 0.00 |
| Starcode MP | 532.00 | 560.00 | 99.45408 | 193.39 | 99.981789 | 4552.00 | 0.00 |
| Bartender | 496.00 | 795.00 | 99.35546 | 41.77 | 99.935342 | 16162.00 | 0.00 |

Simulations: means over all 60 libraries; time is median. Milo: one observation. Shepherd was timed in a later session.

## Supplementary paired accuracy

| Design | LV − Shepherd F1, percentage points | t lower bound | Bootstrap lower | Both > 0 |
|---|---:|---:|---:|---|
| random_mixed | +0.016172 | +0.011819 | +0.012066 | Yes |
| anchored_mixed | +0.014868 | +0.009718 | +0.009941 | Yes |

One-sided alpha 0.05/8 lower bounds; 100,000 whole-library bootstrap resamples. These intervals are supplementary because the Shepherd configuration was selected after the original benchmark. They do not establish universal superiority or retroactively complete the original confirmatory contrasts.

Shepherd completed 121/121 calls. Original +4/automatic results, including all ten failures, remain in the original publication worktree. No final test seed was used for Barbac tuning.

## Development recall

| Input | Total FN | Zero-read truth | True sequence unobserved despite positive reads | Observed identity merged/unassigned |
|---|---:|---:|---:|---:|
| random_substitutions | 124 | 102 | 8 | 14 |
| random_mixed | 125 | 102 | 13 | 10 |
| anchored_substitutions | 119 | 102 | 8 | 9 |
| anchored_mixed | 125 | 102 | 8 | 15 |
| milos | 471 | 409 | 30 | 32 |

These are development references, not final-test decompositions. Most positive misses are one-read errors or tied pairs. More permissive merging cannot identify sequences with no evidence and can erase genuine nearby identities. The previous rejected paired-indel candidate remains excluded.
