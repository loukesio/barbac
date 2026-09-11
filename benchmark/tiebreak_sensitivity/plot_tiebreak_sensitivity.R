#!/usr/bin/env Rscript
# Figure: how large is the between-method difference compared with the spread
# barbac shows across equally valid orderings of its count ties?
#
# Each panel shows one condition and one metric. The blue dots are barbac run
# under 15 alternative tie orders (tie_break = "hash"); the orange diamond is
# the shipped default (tie_break = "sequence"). Competing methods are drawn as
# labelled points on their own rows. Where a competitor falls outside the panel
# its value is
# printed at the margin rather than rescaling the axis, which would flatten the
# spread the figure exists to show.
#
# Usage:
#   Rscript plot_tiebreak_sensitivity.R [out.png]

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                         value = TRUE)[1]))
out_png <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[[1]] else
  file.path(here, "tiebreak_sensitivity.png")

seeds <- readr::read_csv(file.path(here, "tiebreak_seeds.csv"),
                         show_col_types = FALSE)
refs  <- readr::read_csv(file.path(here, "reference_methods.csv"),
                         show_col_types = FALSE)

COND_LABELS <- c(
  sub_only             = "Random 20 bp, no indels",
  low_indel            = "Random 20 bp, 0.5% indels",
  sub_only_structured  = "Anchored 28 bp, no indels",
  low_indel_structured = "Anchored 28 bp, 0.5% indels"
)
METRIC_LABELS <- c(fn_pct = "False negatives (%)",
                   ws_pct = "Wrong-split centroids (%)")

long <- function(d) {
  d %>% tidyr::pivot_longer(c(fn_pct, ws_pct),
                            names_to = "metric", values_to = "value")
}
seeds_l <- long(seeds)
refs_l  <- long(refs)

panel_of <- function(cond, metric) {
  factor(paste0(COND_LABELS[cond], "\n", METRIC_LABELS[metric]),
         levels = as.vector(t(outer(COND_LABELS, METRIC_LABELS,
                                    function(a, b) paste0(a, "\n", b)))))
}
seeds_l$panel <- panel_of(seeds_l$condition, seeds_l$metric)
refs_l$panel  <- panel_of(refs_l$condition,  refs_l$metric)

# Panel limits are set by barbac's own spread, padded. A competitor beyond that
# range is reported as text at the edge: rescaling to fit a 500%-error method
# would compress every barbac point into a single pixel.
lims <- seeds_l %>%
  group_by(panel) %>%
  summarise(lo = min(value), hi = max(value), .groups = "drop") %>%
  mutate(pad = pmax((hi - lo) * 1.6, hi * 0.35, 0.05),
         lo = pmax(0, lo - pad), hi = hi + pad)

# One row per method keeps every label on its own baseline, so nothing can
# collide however close two methods score.
ROWS <- c(barbac = 4, shepherd = 3, starcode = 2, bartender = 1)

refs_l <- refs_l %>%
  left_join(lims, by = "panel") %>%
  mutate(inside = value >= lo & value <= hi,
         y      = ROWS[method],
         label  = sprintf("%.2f%%", value))

inside  <- dplyr::filter(refs_l, inside)
outside <- dplyr::filter(refs_l, !inside) %>%
  mutate(label = sprintf("%.1f%%  →", value))

BLUE   <- "#2a78d6"   # categorical slot 1 - alternative tie orders
ORANGE <- "#eb6834"   # categorical slot 2 - the shipped default
INK    <- "#0b0b0b"
INK2   <- "#52514e"

defaults <- dplyr::filter(seeds_l, is_default)
hashes   <- dplyr::filter(seeds_l, !is_default)

# Span of barbac's tie-break spread, drawn behind the dots as the reference
# interval the competing points are meant to be read against.
spread <- hashes %>%
  group_by(panel) %>%
  summarise(lo = min(value), hi = max(value), .groups = "drop")

p <- ggplot() +
  geom_segment(data = spread,
               aes(x = lo, xend = hi, y = ROWS[["barbac"]], yend = ROWS[["barbac"]]),
               colour = BLUE, linewidth = 2.6, alpha = 0.16,
               lineend = "round") +
  # barbac across alternative tie orders
  geom_point(data = hashes, aes(x = value, y = ROWS[["barbac"]]),
             colour = BLUE, size = 1.8, alpha = 0.8,
             position = position_jitter(height = 0.13, width = 0, seed = 1)) +
  # the shipped default
  geom_point(data = defaults, aes(x = value, y = ROWS[["barbac"]]),
             colour = ORANGE, fill = ORANGE, shape = 23, size = 3.0,
             stroke = 0.5) +
  geom_text(data = defaults, aes(x = value, y = ROWS[["barbac"]] + 0.42,
                                 label = "shipped default"),
            colour = ORANGE, size = 2.4, vjust = 0) +
  # competing methods, one per row
  geom_point(data = inside, aes(x = value, y = y),
             colour = INK, size = 2.0) +
  geom_text(data = inside, aes(x = value, y = y - 0.34, label = label),
            colour = INK2, size = 2.4, vjust = 1) +
  geom_text(data = outside, aes(x = hi, y = y, label = label),
            colour = INK2, size = 2.4, hjust = 1) +
  geom_blank(data = lims, aes(x = lo)) +
  geom_blank(data = lims, aes(x = hi)) +
  facet_wrap(~ panel, ncol = 2, scales = "free_x") +
  scale_y_continuous(
    limits = c(0.45, 5.0),
    breaks = unname(ROWS),
    labels = c("barbac", "Shepherd", "Starcode", "Bartender"),
    name   = NULL
  ) +
  labs(
    x = NULL,
    title = "Between-method differences are comparable to barbac's own tie-break spread",
    subtitle = paste("Blue: barbac under 15 alternative orderings of count-tied barcodes;",
                     "orange: the shipped default ordering.\nBlack: competing methods, one",
                     "run each. Identical data and parameters throughout \u2014 within the blue",
                     "band, only the arbitrary tie order differs.\nValues past a panel's range",
                     "are printed at its right edge."),
    caption = "barbac super_cluster2(distance = 3, merge_ratio = 20), Levenshtein mode. Wrong-split counts false-positive centroids within distance 3 of a true barcode."
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "#f0f0ec", linewidth = 0.25),
    panel.grid.major.x = element_line(colour = "#e8e8e4", linewidth = 0.25),
    strip.text         = element_text(size = 7.6, colour = INK, hjust = 0,
                                      lineheight = 1.25,
                                      margin = margin(b = 3, t = 5)),
    plot.title         = element_text(size = 11, colour = INK, face = "bold",
                                      margin = margin(b = 4)),
    plot.subtitle      = element_text(size = 8, colour = INK2, lineheight = 1.3,
                                      margin = margin(b = 9)),
    plot.caption       = element_text(size = 6.6, colour = INK2, hjust = 0,
                                      margin = margin(t = 8)),
    axis.text.x        = element_text(size = 7, colour = INK2),
    axis.text.y        = element_text(size = 7, colour = INK2, hjust = 1),
    panel.spacing.x    = unit(11, "pt"),
    panel.spacing.y    = unit(9, "pt"),
    plot.margin        = margin(11, 13, 9, 11)
  )

ggsave(out_png, p, width = 8.6, height = 8.6, dpi = 200, bg = "#fcfcfb")
message("wrote ", out_png)

# Numeric summary that the figure is a picture of.
seeds %>%
  group_by(condition) %>%
  summarise(
    default_fn = fn[is_default][1],
    hash_min   = min(fn[!is_default]),
    hash_med   = stats::median(fn[!is_default]),
    hash_max   = max(fn[!is_default]),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
