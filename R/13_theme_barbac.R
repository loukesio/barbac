#' A publication-ready ggplot2 theme for barbac plots
#'
#' A minimal theme with a solid inner border, visible tick marks on both
#' axes, and a centred bold title. Designed to complement
#' [barbac_ts_area()] and other lineage-visualisation output.
#'
#' Font handling is deliberately conservative. The default \code{family}
#' is \code{"sans"} so the theme works on every platform without extra
#' setup. Pass \code{family = "Avenir"} if you know the font is installed
#' locally (macOS ships with it); pass one of the Google Fonts that the
#' \pkg{showtext} package can register at runtime (e.g. \code{"Nunito
#' Sans"}, \code{"Inter"}, \code{"Montserrat"}) if you want a similar
#' geometric feel without a local font install. See
#' \code{vignette("barbac")} for an example.
#'
#' The commercial Avenir font is licensed by Adobe/Linotype and cannot
#' be redistributed inside an R package; the theme therefore never
#' bundles it.
#'
#' @param base_size Base font size in points. Default 14.
#' @param family Font family. Default \code{"sans"}. Pass
#'   \code{"Avenir"} for the macOS default, or any Google Font name
#'   after registering it with \code{sysfonts::font_add_google()}.
#' @param border_colour Colour of the panel border and axis ticks.
#'   Default \code{"#333333"}.
#'
#' @return A ggplot2 theme object; add it to any plot with \code{+}.
#'
#' @examples
#' \dontrun{
#' # Straight-forward use with sans (works everywhere).
#' p <- barbac_ts_area(my_time_series_df)
#' p + theme_barbac()
#'
#' # Using macOS Avenir if you have it.
#' p + theme_barbac(family = "Avenir")
#'
#' # Using a Google Font via showtext.
#' sysfonts::font_add_google("Nunito Sans", "nunito")
#' showtext::showtext_auto()
#' p + theme_barbac(family = "nunito")
#' }
#'
#' @export
#' @importFrom ggplot2 theme_minimal theme element_text element_line
#'   element_rect element_blank unit
theme_barbac <- function(base_size    = 14,
                         family       = "sans",
                         border_colour = "#333333") {
  ggplot2::theme_minimal(base_size = base_size, base_family = family) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5,
                                               size  = base_size * 1.2,
                                               face  = "bold",
                                               family = family),
      legend.position  = "bottom",
      axis.ticks.x     = ggplot2::element_line(colour = border_colour),
      axis.ticks.y     = ggplot2::element_line(colour = border_colour),
      axis.ticks.length = ggplot2::unit(0.26, "cm"),
      panel.border     = ggplot2::element_rect(colour = border_colour, fill = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      text             = ggplot2::element_text(family = family),
      axis.text        = ggplot2::element_text(size = base_size * 1.25),
      axis.title       = ggplot2::element_text(size = base_size * 1.45)
    )
}
