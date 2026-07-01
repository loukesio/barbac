# Run Full Barbac CLI Pipeline

This function runs the full CLI pipeline for barcode analysis: FastQC
-\> PEAR -\> Minimap2 -\> BAM Stats.

## Usage

``` r
run_cli_pipeline(
  sample_table,
  reference,
  output_dir = "results",
  verbose = TRUE,
  log_file = NULL,
  create_output_dir = TRUE
)
```

## Arguments

- sample_table:

  Path to samples.csv or a data.frame/tibble containing \`sample\`,
  \`R1\`, and optionally \`R2\`.

- reference:

  Path to reference FASTA file.

- output_dir:

  Directory to write output files. Default: "results".

- verbose:

  Whether to print progress messages to console. Default: TRUE.

- log_file:

  File path for logging. If NULL, uses output_dir/pipeline.log. Default:
  NULL.

- create_output_dir:

  Whether to create output_dir if it doesn't exist. Default: TRUE.

## Value

Invisible list containing:

- commands: Character vector of system commands executed

- output_dir: Path to output directory

- stats: Data frame with BAM statistics

## Examples

``` r
if (FALSE) { # \dontrun{
# Run with default output directory (./results)
run_cli_pipeline(
  sample_table = "samples.csv",
  reference = "barcodes.fasta"
)

# Run with custom output directory
run_cli_pipeline(
  sample_table = "samples.csv",
  reference = "barcodes.fasta",
  output_dir = "/path/to/my_results"
)
} # }
```
