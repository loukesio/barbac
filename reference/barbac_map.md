# Map reads to a reference (preferred name for \`run_cli_pipeline()\`)

\`barbac_map()\` forwards all arguments to \[run_cli_pipeline()\], which
runs FastQC, optional PEAR merging, minimap2 and samtools and writes
sorted, indexed BAM files with mapping statistics. Both names are
supported.

## Usage

``` r
barbac_map(...)
```

## Arguments

- ...:

  Arguments passed to \[run_cli_pipeline()\], such as the sample table,
  reference and output directory.

## Value

The invisible list returned by \[run_cli_pipeline()\].

## See also

\[run_cli_pipeline()\] for the full argument documentation.
