# Stacked-area plot of barcode lineage frequencies over time

Prepares a Bartender-wide or long-format barcode count table into
per-timepoint frequencies and plots a stacked-area chart of lineage
trajectories over time. Barcodes that first appear at later timepoints
can optionally be carried through the full series; missing (barcode,
timepoint) cells are filled either with a small epsilon or with zero.

## Usage

``` r
barbac_ts_area(
  data,
  id_col = "barcode",
  time_col = "time",
  count_col = "counts",
  wide_time_prefix = "time_point_",
  min_total_count = 10,
  include_late = TRUE,
  fill_missing = c("epsilon", "zero"),
  epsilon = 1e-06,
  palette = NULL,
  time_zero_shift = FALSE,
  x_breaks = NULL,
  title = NULL,
  x_lab = "Time",
  y_lab = "Barcode frequency",
  show_legend = FALSE,
  theme = theme_barbac(),
  interactive = FALSE
)
```

## Arguments

- data:

  A `data.frame`/`tibble` or a file path to a CSV.

- id_col:

  Character. Name of the lineage identifier column in long-format input
  (default `"barcode"`); used as the destination name when a
  Bartender-wide `Cluster.ID` column is renamed.

- time_col:

  Character. Name of the time column (default `"time"`).

- count_col:

  Character. Name of the count column (default `"counts"`).

- wide_time_prefix:

  Character prefix used to detect and parse Bartender-style wide time
  columns (default `"time_point_"`).

- min_total_count:

  Numeric. Lineages whose summed count over all timepoints is below this
  threshold are dropped (default `10`).

- include_late:

  Logical. If `TRUE` (default), keep lineages that first appear after
  the earliest timepoint and complete missing earlier timepoints with
  `fill_missing`. If `FALSE`, only lineages present at the earliest
  timepoint are plotted.

- fill_missing:

  How to fill missing (barcode, timepoint) cells after grid completion.
  One of `"epsilon"` (default) or `"zero"`.

- epsilon:

  Numeric. Value used when `fill_missing = "epsilon"` (default `1e-6`).

- palette:

  Optional character vector of colours. If `NULL` (default), a
  continuous `PNWColors::pnw_palette("Sailboat")` is expanded to the
  number of lineages.

- time_zero_shift:

  Logical. If `TRUE`, remap `time == 1` to `0`. Only useful when the
  input labels the earliest timepoint as 1 (e.g. Bartender
  `time_point_1` columns) and you want the plot to start at 0. Default
  `FALSE`. Passing `TRUE` on data that already contains an explicit
  `time == 0` row will collapse those rows with the original `time == 1`
  rows, and counts within each (barcode, shifted-time) cell are summed.

- x_breaks:

  Optional numeric vector of x-axis breaks. If `NULL`, uses integer
  steps across the time range.

- title, x_lab, y_lab:

  Plot text (defaults sensible for barcode frequency plots).

- show_legend:

  Logical. If `TRUE`, show the barcode legend. With many lineages this
  is usually too dense to be useful, so defaults to `FALSE`.

- theme:

  A ggplot2 theme object added to the plot before it is returned.
  Default
  [`theme_barbac()`](https://loukesio.github.io/barbac/reference/theme_barbac.md).
  Pass `NULL` for an unthemed plot, or any other `ggplot2::theme_*` to
  swap.

- interactive:

  One of `FALSE` (static, default), `TRUE` / `"ggiraph"` (recommended
  interactive backend — preserves the exact ggplot look and produces a
  much lighter HTML than plotly for many-band plots), or `"plotly"`
  (uses
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html) —
  offers built-in zoom / pan / lasso but translation drifts slightly
  from the ggplot theme). The interactive modes require the
  corresponding package to be installed.

## Value

When `interactive = FALSE`, a `ggplot` object. When `interactive = TRUE`
or `"ggiraph"`, a `ggiraph` htmlwidget. When `"plotly"`, a `plotly`
htmlwidget.

## Details

Two input layouts are auto-detected:

- **Long**: one row per (barcode, timepoint) pair with columns `id_col`,
  `time_col`, and `count_col`.

- **Bartender wide**: one row per lineage with an id column `Cluster.ID`
  and abundance columns whose names begin with `wide_time_prefix`
  (default `"time_point_"`). Any `Center` and `Cluster.Score` columns
  are dropped.

If neither schema matches, the function stops with a message that lists
the expected column names.

## Examples

``` r
if (FALSE) { # \dontrun{
# long-format input
df <- tibble::tibble(
  barcode = rep(c("A", "B", "C"), each = 5),
  time    = rep(0:4, times = 3),
  counts  = c(100, 80, 60, 40, 20,
               20, 40, 60, 80, 100,
                0,  0, 10, 50, 40)
)
barbac_ts_area(df, min_total_count = 0)

# Bartender wide-format input auto-detects
bt <- tibble::tibble(
  Cluster.ID    = c("bc1", "bc2"),
  Center        = c("ACGT", "ACGA"),
  Cluster.Score = c(0.9, 0.8),
  time_point_1  = c(100, 20),
  time_point_2  = c(80,  40),
  time_point_3  = c(60,  60)
)
barbac_ts_area(bt, min_total_count = 0)
} # }
```
