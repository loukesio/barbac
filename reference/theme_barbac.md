# A publication-ready ggplot2 theme for barbac plots

A minimal theme with a solid inner border, visible tick marks on both
axes, and a centred bold title. Designed to complement
\[barbac_ts_area()\] and other lineage-visualisation output.

## Usage

``` r
theme_barbac(base_size = 14, family = "sans", border_colour = "#333333")
```

## Arguments

- base_size:

  Base font size in points. Default 14.

- family:

  Font family. Default `"sans"`. Pass `"Avenir"` for the macOS default,
  or any Google Font name after registering it with
  [`sysfonts::font_add_google()`](https://rdrr.io/pkg/sysfonts/man/font_add_google.html).

- border_colour:

  Colour of the panel border and axis ticks. Default `"#333333"`.

## Value

A ggplot2 theme object; add it to any plot with `+`.

## Details

Font handling is deliberately conservative. The default `family` is
`"sans"` so the theme works on every platform without extra setup. Pass
`family = "Avenir"` if you know the font is installed locally (macOS
ships with it); pass one of the Google Fonts that the showtext package
can register at runtime (e.g. `"Nunito Sans"`, `"Inter"`,
`"Montserrat"`) if you want a similar geometric feel without a local
font install. See
[`vignette("barbac")`](https://loukesio.github.io/barbac/articles/barbac.md)
for an example.

The commercial Avenir font is licensed by Adobe/Linotype and cannot be
redistributed inside an R package; the theme therefore never bundles it.

## Examples

``` r
if (FALSE) { # \dontrun{
# Straight-forward use with sans (works everywhere).
p <- barbac_ts_area(my_time_series_df)
p + theme_barbac()

# Using macOS Avenir if you have it.
p + theme_barbac(family = "Avenir")

# Using a Google Font via showtext.
sysfonts::font_add_google("Nunito Sans", "nunito")
showtext::showtext_auto()
p + theme_barbac(family = "nunito")
} # }
```
