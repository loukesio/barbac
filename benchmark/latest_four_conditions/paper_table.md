**Table 2. Latest barcode clustering performance at distance three.** R-S: random barcodes with substitutions; R-I: random barcodes with substitutions and indels; A-S and A-I: the corresponding anchored designs. Milos denotes the Johnson et al. (2023) reference simulation. FN and FP are absolute barcode counts; for R-S, R-I, A-S and A-I, counts and centroid F1 are means across seeds 42, 43 and 44 (10,000 true barcodes and one million reads per seed). The Milos row for each method uses one fixed dataset of 100,000 true barcodes and 24,996,128 reads. Bold F1 values identify the highest score within a dataset. F1 measures exact centroid recovery and differs from the cluster-label F1 used in Section 3.1. Time is one serial workflow observation in seconds on the same development machine: seed 42 for the four smaller simulations and the full dataset for Milos. It includes program startup, required format conversion, clustering, and centroid/member exports; shared input staging and scoring are excluded. These timings have no confidence intervals. Both barbac modes use native v13 and support ordering; LV additionally enables the experimental Poisson indel option. Starcode sphere and default message passing (MP) are reported separately. All methods use distance three, with its method-specific Hamming or Levenshtein interpretation.

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
