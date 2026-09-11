# Generate Summary Statistics and Plots for Extracted Barcodes

Generate Summary Statistics and Plots for Extracted Barcodes

## Usage

``` r
barbac_xtr.stats(
  file,
  barcode_length,
  fill_color = "#FF5349",
  save_plot = FALSE,
  plot_width = 12,
  plot_height = 8,
  verbose = TRUE,
  panel_labels = FALSE,
  return_details = FALSE
)
```

## Arguments

- file:

  A data frame or path to CSV file with columns: barcode, counts,
  barcode_length

- barcode_length:

  A numeric vector of length 2 specifying desired barcode length range
  \[min, max\]

- fill_color:

  Histogram fill color. Default is "#FF5349" (red).

- save_plot:

  Logical or character. If TRUE, saves to default name. If character,
  saves to that path. Default is FALSE.

- plot_width:

  Numeric. Width of saved plot in inches. Default is 12.

- plot_height:

  Numeric. Height of saved plot in inches. Default is 8.

- verbose:

  Logical. Print summary statistics. Default is TRUE.

- panel_labels:

  Logical. Add A, B, C, D labels to the combined diagnostic.

- return_details:

  Logical. Return the combined plot, individual histogram plots, numeric
  length summary and per-sequence entropy for report composition. The
  default preserves the patchwork return value.

## Value

A patchwork plot, or with \`return_details = TRUE\`, a list containing
\`plot\`, \`plots\`, \`length_summary\` and \`sequence_entropy\`. The
numeric summary can be passed to \`gt::gt()\` for an HTML report.
Entropy is measured in bits across nucleotide frequencies within each
sequence, not across lineages.
