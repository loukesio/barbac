# Report which external command-line tools barbac can find

The CLI pipeline shells out to FastQC, MultiQC, PEAR, minimap2 and
samtools. They may come from the conda environment
[`configure_environment()`](https://loukesio.github.io/barbac/reference/configure_environment.md)
creates, or already be installed system-wide; either works, so both are
searched and the result says which one supplied each tool.

## Usage

``` r
check_barbac_tools(
  env_name = "barbac_env",
  tools = c("fastqc", "multiqc", "pear", "minimap2", "samtools")
)
```

## Arguments

- env_name:

  Conda environment to search first. Default: "barbac_env".

- tools:

  Character vector of tool names to look for.

## Value

A data.frame with one row per tool: `tool`, `available`, `source` (the
environment name, `"PATH"`, or NA when missing), `version` (NA when the
tool reports none), and `path`.
