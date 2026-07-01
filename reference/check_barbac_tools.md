# Check availability of tools in barbac conda environment

Check availability of tools in barbac conda environment

## Usage

``` r
check_barbac_tools(
  env_name = "barbac_env",
  tools = c("fastqc", "multiqc", "pear", "minimap2", "samtools")
)
```

## Arguments

- env_name:

  Name of the conda environment (default: "barbac_env")

- tools:

  Character vector of tool names to check

## Value

A data.frame with tool names and availability
