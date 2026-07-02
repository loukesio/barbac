# Per-sample bar plot of mapping counts

Companion visualisation for \[summarise_bam_stats()\]. Draws a bar plot
of mapped vs unmapped read counts across samples so outlier samples are
visible at a glance.

## Usage

``` r
plot_bam_stats(
  stats,
  position = c("stack", "dodge"),
  mapped_color = "#2C4A63",
  unmapped_color = "#B0B7BF"
)
```

## Arguments

- stats:

  Either a data.frame with columns \`sample\`, \`mapped\`, \`unmapped\`
  (as returned by \[summarise_bam_stats()\]) or a path to a directory of
  sorted BAM files (in which case \`summarise_bam_stats()\` is called
  first).

- position:

  Character. Bar layout: \`"stack"\` (default) or \`"dodge"\`.

- mapped_color:

  Character. Fill colour for mapped reads. Default \`"#2C4A63"\` (barbac
  primary).

- unmapped_color:

  Character. Fill colour for unmapped reads. Default \`"#B0B7BF"\`.

## Value

A \`ggplot\` object.
