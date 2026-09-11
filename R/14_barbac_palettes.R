# LTC palettes by Loukas Theodosiou, using the curated fills from ggvmap.
# See inst/COPYRIGHTS for source attribution and the MIT permission notice.
.barbac_ltc_palettes <- list(
  paloma     = c("#83AF9B", "#C8C8A9", "#f8da8a", "#f7bf95", "#fe8ca1"),
  maya       = c("#3d5a80", "#98c1d9", "#ee6c4d", "#293241"),
  dora       = c("#52777A", "#542437", "#C02942", "#D95B43", "#ECD078"),
  ploen      = c("#3F5671", "#83A1C3", "#CEB5C8", "#FAC898", "#B17776"),
  olga       = c("#c9e3c2", "#8bc8cb", "#eccd80", "#f5ab70", "#9c87a1"),
  mterese    = c("#f7ddaa", "#fac3ad", "#f897a1", "#9298BA", "#9cbeed"),
  gaby       = c("#fceaab", "#f1a890", "#a8c4cc", "#82A0C2", "#85496F"),
  franscoise = c("#5980B1", "#b96a8d", "#A55062", "#E05256", "#E9A986"),
  fernande   = c("#ff7676", "#F9D662", "#7cab7d", "#75B7D1"),
  sylvie     = c("#E8B961", "#E88170", "#C6BDE8", "#5DB7C4", "#FD95BC"),
  expevo     = c("#FC4E07", "#E7B800", "#00AFBB", "#8B4769", "#1d457f", "#808080"),
  minou      = c("#00798c", "#d1495b", "#edae49", "#66a182", "#2e4057", "#8d96a3"),
  kiss       = c("#FF7C7E", "#FEC300", "#9E3F71", "#31BCBA", "#E20035"),
  hat        = c("#efb306", "#eb990c", "#e8351e", "#cd023d", "#852f88",
                 "#4e54ac", "#0f8096", "#7db954", "#17a769"),
  reading    = c("#EFBC68", "#919F89", "#EDBDAE", "#57717C", "#5F97A4",
                 "#CAEAC8", "#95A1AE", "#C8CFD6"),
  alger      = c("#1A5B5B", "#ACC8BE", "#F4AB5C", "#D1422F"),
  trio1      = c("#0E7175", "#FD7901", "#C35BCA"),
  trio2      = c("#89973D", "#E8B92F", "#A45E41"),
  trio3      = c("#E69F00", "#56B4E9", "#009E73"),
  trio4      = c("#94475E", "#364C54", "#E5A11F"),
  heatmap0   = c("#001219", "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
                 "#EE9B00", "#CA6702", "#AE2012", "#9B2226"),
  pantone23  = c("#7A92A5", "#1F2C43", "#FFB000", "#842c48", "#46483d"),
  remains    = c("#69326E", "#FF6D1F", "#EED455"),
  midnight   = c("#16232A", "#FF5B04", "#075056"),
  lincoln    = c("#C9C1B1", "#2C3B4D", "#FFB162", "#A35139", "#1B2632"),
  luminaries = c("#FF5B04", "#075056", "#233038", "#F4D47C", "#D3DBDD"),
  seafarer   = c("#013D5A", "#BDD3CE", "#708C69", "#E4A25B"),
  shuggie    = c("#5B5F8D", "#9BB29E", "#DA6B51", "#F1DCBA", "#484149"),
  heatmap1   = c("#4d7799", "#7fa4c4", "#c5c8d4", "#d48e95", "#b5515b"),
  heatmap2   = c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0"),
  heatmap3   = c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6"),
  casa_natal = c("#245E55", "#ED773C", "#808BC5", "#C63F3E", "#EAC119",
                 "#EAA7C7", "#9ED6DF")
)

#' Built-in LTC colour palettes
#'
#' The 32 LTC palettes by Loukas Theodosiou are included directly in barbac;
#' no additional palette package is needed. These use the curated fills from
#' ggvmap, which omit selected black and near-white entries from \code{alger},
#' \code{hat}, \code{casa_natal}, \code{luminaries}, \code{seafarer},
#' \code{lincoln}, \code{midnight}, \code{remains}, and \code{maya}.
#'
#' \code{barbac_palettes()} returns all base colour vectors.
#' \code{barbac_palette()} selects one palette and optionally interpolates
#' it end to end. Names ignore case, spaces, underscores and hyphens:
#' \code{"Casa Natal"} and \code{"casa_natal"} are equivalent.
#' Interpolated colours identify lineages, not their abundance. With many
#' lineages, neighbouring colours can be visually indistinguishable.
#'
#' @param name A built-in LTC palette name, such as \code{"alger"}.
#' @param n Number of colours to generate by RGB interpolation, including
#'   both endpoints when \code{n >= 2}. \code{NULL} returns the original
#'   base colours; zero returns an empty vector and one returns the first colour.
#' @return \code{barbac_palettes()} returns a named list of colour vectors;
#'   \code{barbac_palette()} returns a character vector of colours.
#' @seealso \code{\link{barbac_ts_area}}
#' @examples
#' names(barbac_palettes())
#' barbac_palette("alger")
#' barbac_palette("Casa Natal", n = 100)
#' @export
barbac_palettes <- function() {
  .barbac_ltc_palettes
}

.barbac_palette_key <- function(x) gsub("[ _-]", "", tolower(x))

#' @rdname barbac_palettes
#' @export
barbac_palette <- function(name, n = NULL) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("`name` must be a single LTC palette name. See names(barbac_palettes()).",
         call. = FALSE)
  }
  index <- match(.barbac_palette_key(name),
                 .barbac_palette_key(names(.barbac_ltc_palettes)))
  if (is.na(index)) {
    stop("Unknown LTC palette: ", name,
         ". See names(barbac_palettes()) for available names.", call. = FALSE)
  }
  colours <- .barbac_ltc_palettes[[index]]
  if (is.null(n)) return(colours)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 0 || n != floor(n) || n > .Machine$integer.max) {
    stop("`n` must be a single non-negative whole number.", call. = FALSE)
  }
  if (n == 0) return(character())
  grDevices::colorRampPalette(colours)(n)
}

# Internal: retain the default and custom-vector behaviour of barbac_ts_area.
.barbac_ts_colours <- function(palette, n) {
  if (is.null(palette)) {
    base <- PNWColors::pnw_palette("Sailboat", n = n, type = "continuous")
    return(grDevices::colorRampPalette(base)(n))
  }
  if (!is.character(palette) || !length(palette) || anyNA(palette)) {
    stop("`palette` must be an LTC palette name or a non-empty colour vector.",
         call. = FALSE)
  }
  if (length(palette) == 1L && .barbac_palette_key(palette) %in%
      .barbac_palette_key(names(.barbac_ltc_palettes))) {
    return(barbac_palette(palette, n))
  }
  tryCatch(grDevices::col2rgb(palette), error = function(e) {
    stop("Unknown palette name or invalid colour in `palette`. ",
         "See names(barbac_palettes()) or supply valid R colours.", call. = FALSE)
  })
  if (length(palette) < n) {
    grDevices::colorRampPalette(palette)(n)
  } else {
    palette[seq_len(n)]
  }
}
