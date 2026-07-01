# Fast Centroid-Based Sequence Clustering

Abundance-ranked centroid clustering with likelihood-based best-parent
selection and a distance-aware count-ratio merge guard. Supports both
Levenshtein (handles indels) and Hamming distances.

## Usage

``` r
super_cluster2(
  input_path,
  distance = 3,
  method = c("lv", "hamming", "osa", "dl", "lcs", "qgram", "cosine", "jaccard", "jw",
    "soundex"),
  barcode_col = "barcode",
  counts_col = "counts",
  output_dir = NULL,
  file_pattern = "\\.csv$",
  verbose = TRUE,
  use_cpp = TRUE,
  use_kmer_filter = TRUE,
  kmer_size = 5L,
  min_shared_kmers = 2L,
  merge_ratio = 20,
  error_rate = 0.005,
  q = 2
)
```

## Arguments

- input_path:

  Character string or data.frame.

- distance:

  Numeric. Maximum edit distance. Default: 3.

- method:

  Character string. "lv" (Levenshtein, default) or "hamming".

- barcode_col:

  Character string. Default: "barcode".

- counts_col:

  Character string. Default: "counts".

- output_dir:

  Character string or NULL. Default: NULL.

- file_pattern:

  Character string. Default: "\\csv\$".

- verbose:

  Logical. Default: TRUE.

- use_cpp:

  Logical. Default: TRUE.

- use_kmer_filter:

  Logical. Default: TRUE.

- kmer_size:

  Integer. Seed size for LV index. Default: 5.

- min_shared_kmers:

  Integer. Kept for API compatibility. Default: 2.

- merge_ratio:

  Numeric. Base count-ratio for the distance-aware merge guard.
  Effective ratio increases with distance. Default: 20.

- error_rate:

  Numeric. Approximate per-base error rate for likelihood scoring.
  Default: 0.005.

- q:

  Integer. Q-gram size for qgram/jaccard/cosine. Default: 2.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
columns: cluster_id, central_barcode, all_barcodes, all_counts,
sum_counts.
