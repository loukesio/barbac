# Configure Conda Environment for barbac

Configure Conda Environment for barbac

## Usage

``` r
configure_environment(
  env_name = "barbac_env",
  channels = c("conda-forge", "bioconda"),
  confirm = TRUE
)
```

## Arguments

- env_name:

  Name of the conda environment to create. Default is "barbac_env".

- channels:

  List of conda channels to use (default: c("conda-forge", "bioconda"))

- confirm:

  If FALSE, it will skip the interactive menu. Default is TRUE.

## Value

Prints progress messages and returns the environment path if successful.
