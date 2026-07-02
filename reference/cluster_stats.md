# Post-clustering QC summary

Reports high-level diagnostics from \[super_cluster2()\] output to help
sanity-check \`distance\` and \`merge_ratio\` choices before running
downstream time-series analysis. A very high singleton fraction usually
means \`distance\` was too tight; a dominant single cluster with a small
\`top1_frac\` \<\< expected can indicate over-aggressive merging.

## Usage

``` r
cluster_stats(
  clusters,
  top_n = c(1, 10, 100),
  singleton_max = 1L,
  verbose = TRUE
)
```

## Arguments

- clusters:

  Tibble returned by \[super_cluster2()\], or a path to a CSV of that
  shape. Must contain a \`sum_counts\` column.

- top_n:

  Integer vector. Rank thresholds for cumulative-abundance share.
  Default \`c(1, 10, 100)\`.

- singleton_max:

  Integer. Maximum \`sum_counts\` treated as a singleton. Default
  \`1L\`.

- verbose:

  Logical. Print the summary to console. Default \`TRUE\`.

## Value

A one-row tibble with columns:

- \`n_clusters\`

- \`total_reads\`

- \`n_singletons\`, \`singleton_frac\`

- \`largest_cluster_size\`, \`largest_cluster_frac\`

- \`topK_frac\` for each \`K\` in \`top_n\`

## Examples

``` r
if (FALSE) { # \dontrun{
result <- super_cluster2("sample1_barcodes.csv", distance = 3)
cluster_stats(result)
} # }
```
