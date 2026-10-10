# Cluster barcodes (preferred name for \`super_cluster2()\`)

\`barbac_cluster()\` is the public name for the native abundance-aware
centroid clustering engine. It forwards all arguments to
\[super_cluster2()\] and returns the same object. Both names are
supported.

## Usage

``` r
barbac_cluster(...)
```

## Arguments

- ...:

  Arguments passed to \[super_cluster2()\], such as \`input_path\`,
  \`distance\`, \`method\`, \`merge_ratio\`, \`error_rate\`,
  \`tie_break\` and \`indel_model\`.

## Value

The tibble returned by \[super_cluster2()\]: one row per cluster with
\`cluster_id\`, \`central_barcode\`, \`all_barcodes\`, \`all_counts\`
and \`sum_counts\`.

## See also

\[super_cluster2()\] for the full argument documentation,
\[cluster_stats()\] for post-clustering diagnostics.

## Examples

``` r
input <- data.frame(
  barcode = c("AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAT"),
  counts = c(1000, 10)
)
barbac_cluster(input, distance = 3, verbose = FALSE)
#> # A tibble: 1 × 5
#>   cluster_id central_barcode      all_barcodes all_counts sum_counts
#>   <chr>      <chr>                <list>       <list>          <int>
#> 1 group1     AAAAAAAAAAAAAAAAAAAA <chr [2]>    <int [2]>        1010
```
