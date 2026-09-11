# Built-in LTC colour palettes

The 32 LTC palettes by Loukas Theodosiou are included directly in
barbac; no additional palette package is needed. These use the curated
fills from ggvmap, which omit selected black and near-white entries from
`alger`, `hat`, `casa_natal`, `luminaries`, `seafarer`, `lincoln`,
`midnight`, `remains`, and `maya`.

## Usage

``` r
barbac_palettes()

barbac_palette(name, n = NULL)
```

## Arguments

- name:

  A built-in LTC palette name, such as `"alger"`.

- n:

  Number of colours to generate by RGB interpolation, including both
  endpoints when `n >= 2`. `NULL` returns the original base colours;
  zero returns an empty vector and one returns the first colour.

## Value

`barbac_palettes()` returns a named list of colour vectors;
`barbac_palette()` returns a character vector of colours.

## Details

`barbac_palettes()` returns all base colour vectors. `barbac_palette()`
selects one palette and optionally interpolates it end to end. Names
ignore case, spaces, underscores and hyphens: `"Casa Natal"` and
`"casa_natal"` are equivalent. Interpolated colours identify lineages,
not their abundance. With many lineages, neighbouring colours can be
visually indistinguishable.

## See also

[`barbac_ts_area`](https://loukesio.github.io/barbac/reference/barbac_ts_area.md)

## Examples

``` r
names(barbac_palettes())
#>  [1] "paloma"     "maya"       "dora"       "ploen"      "olga"      
#>  [6] "mterese"    "gaby"       "franscoise" "fernande"   "sylvie"    
#> [11] "expevo"     "minou"      "kiss"       "hat"        "reading"   
#> [16] "alger"      "trio1"      "trio2"      "trio3"      "trio4"     
#> [21] "heatmap0"   "pantone23"  "remains"    "midnight"   "lincoln"   
#> [26] "luminaries" "seafarer"   "shuggie"    "heatmap1"   "heatmap2"  
#> [31] "heatmap3"   "casa_natal"
barbac_palette("alger")
#> [1] "#1A5B5B" "#ACC8BE" "#F4AB5C" "#D1422F"
barbac_palette("Casa Natal", n = 100)
#>   [1] "#245E55" "#305F53" "#3C6151" "#486250" "#54644E" "#60654D" "#6D674B"
#>   [8] "#79684A" "#856A48" "#916B47" "#9D6D45" "#AA6E44" "#B67042" "#C27141"
#>  [15] "#CE733F" "#DA743E" "#E6763C" "#E97740" "#E37848" "#DC7A50" "#D57B59"
#>  [22] "#CF7C61" "#C87D69" "#C27E71" "#BB807A" "#B48182" "#AE828A" "#A78393"
#>  [29] "#A1849B" "#9A86A3" "#9387AC" "#8D88B4" "#8689BC" "#808AC4" "#8486BC"
#>  [36] "#8881B4" "#8C7DAC" "#9078A4" "#95739C" "#996F93" "#9D6A8B" "#A16683"
#>  [43] "#A6617B" "#AA5C73" "#AE586A" "#B25362" "#B74F5A" "#BB4A52" "#BF454A"
#>  [50] "#C34142" "#C7423C" "#C94A3A" "#CB5238" "#CD5A36" "#CF6233" "#D26A31"
#>  [57] "#D4722F" "#D67A2D" "#D8812A" "#DA8928" "#DC9126" "#DF9924" "#E1A121"
#>  [64] "#E3A91F" "#E5B11D" "#E7B91B" "#EAC019" "#EABF23" "#EABD2E" "#EABC38"
#>  [71] "#EABA43" "#EAB94D" "#EAB758" "#EAB562" "#EAB46D" "#EAB277" "#EAB182"
#>  [78] "#EAAF8D" "#EAAE97" "#EAACA2" "#EAAAAC" "#EAA9B7" "#EAA7C1" "#E7A8C7"
#>  [85] "#E3ABC9" "#DEAECA" "#D9B0CC" "#D5B3CD" "#D0B6CF" "#CCB9D0" "#C7BCD1"
#>  [92] "#C2BFD3" "#BEC2D4" "#B9C4D6" "#B5C7D7" "#B0CAD9" "#ABCDDA" "#A7D0DC"
#>  [99] "#A2D3DD" "#9ED6DF"
```
