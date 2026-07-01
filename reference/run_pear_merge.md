# Merge paired-end reads using PEAR

This function reads a \`samples.csv\` file from the given path and
merges reads using PEAR. If only \`R1\` is available, the sample is
skipped. Output is saved to \`results/merged/\` inside the same folder.

## Usage

``` r
run_pear_merge(sample_csv)
```

## Arguments

- sample_csv:

  Full path to \`samples.csv\`

## Value

Character vector of PEAR commands run.
