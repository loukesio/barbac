# Use barbac conda environment in the current R session

Sets PATH to include the barbac conda environment so CLI tools work in
R.

## Usage

``` r
use_barbac_env(env_name = "barbac_env")
```

## Arguments

- env_name:

  Name of the conda environment. Default is "barbac_env"

## Value

The full path to the environment's bin directory, invisibly
