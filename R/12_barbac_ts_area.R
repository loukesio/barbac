#' Stacked-area plot of barcode lineage frequencies over time
#'
#' Prepares a Bartender-wide or long-format barcode count table into
#' per-timepoint frequencies and plots a stacked-area chart of lineage
#' trajectories over time. Barcodes that first appear at later timepoints
#' can optionally be carried through the full series; missing
#' (barcode, timepoint) cells are filled either with a small epsilon or
#' with zero.
#'
#' Two input layouts are auto-detected:
#' \itemize{
#'   \item \strong{Long}: one row per (barcode, timepoint) pair with
#'     columns \code{id_col}, \code{time_col}, and \code{count_col}.
#'   \item \strong{Bartender wide}: one row per lineage with an id
#'     column \code{Cluster.ID} and abundance columns whose names begin
#'     with \code{wide_time_prefix} (default \code{"time_point_"}).
#'     Any \code{Center} and \code{Cluster.Score} columns are dropped.
#' }
#' If neither schema matches, the function stops with a message that
#' lists the expected column names.
#'
#' @param data A \code{data.frame}/\code{tibble} or a file path to a CSV.
#' @param id_col Character. Name of the lineage identifier column in
#'   long-format input (default \code{"barcode"}); used as the destination
#'   name when a Bartender-wide \code{Cluster.ID} column is renamed.
#' @param time_col Character. Name of the time column (default \code{"time"}).
#' @param count_col Character. Name of the count column
#'   (default \code{"counts"}).
#' @param wide_time_prefix Character prefix used to detect and parse
#'   Bartender-style wide time columns (default \code{"time_point_"}).
#' @param min_total_count Numeric. Lineages whose summed count over all
#'   timepoints is below this threshold are dropped (default \code{10}).
#' @param include_late Logical. If \code{TRUE} (default), keep lineages
#'   that first appear after the earliest timepoint and complete missing
#'   earlier timepoints with \code{fill_missing}. If \code{FALSE}, only
#'   lineages present at the earliest timepoint are plotted.
#' @param fill_missing How to fill missing (barcode, timepoint) cells
#'   after grid completion. One of \code{"epsilon"} (default) or
#'   \code{"zero"}.
#' @param epsilon Numeric. Value used when
#'   \code{fill_missing = "epsilon"} (default \code{1e-6}).
#' @param palette Optional character vector of colours. If \code{NULL}
#'   (default), a continuous \code{PNWColors::pnw_palette("Sailboat")}
#'   is expanded to the number of lineages.
#' @param time_zero_shift Logical. If \code{TRUE}, remap \code{time == 1}
#'   to \code{0}. Only useful when the input labels the earliest
#'   timepoint as 1 (e.g. Bartender \code{time_point_1} columns) and you
#'   want the plot to start at 0. Default \code{FALSE}. Passing
#'   \code{TRUE} on data that already contains an explicit \code{time == 0}
#'   row will collapse those rows with the original \code{time == 1}
#'   rows, and counts within each (barcode, shifted-time) cell are summed.
#' @param x_breaks Optional numeric vector of x-axis breaks. If
#'   \code{NULL}, uses integer steps across the time range.
#' @param title,x_lab,y_lab Plot text (defaults sensible for barcode
#'   frequency plots).
#' @param show_legend Logical. If \code{TRUE}, show the barcode legend.
#'   With many lineages this is usually too dense to be useful, so
#'   defaults to \code{FALSE}.
#' @param theme A ggplot2 theme object added to the plot before it is
#'   returned. Default \code{\link{theme_barbac}()}. Pass \code{NULL}
#'   for an unthemed plot, or any other \code{ggplot2::theme_*} to swap.
#' @param interactive One of \code{FALSE} (static, default),
#'   \code{TRUE} / \code{"ggiraph"} (recommended interactive backend —
#'   preserves the exact ggplot look and produces a much lighter HTML
#'   than plotly for many-band plots), or \code{"plotly"} (uses
#'   \code{plotly::ggplotly()} — offers built-in zoom / pan / lasso but
#'   translation drifts slightly from the ggplot theme). The interactive
#'   modes require the corresponding package to be installed.
#'
#' @return When \code{interactive = FALSE}, a \code{ggplot} object. When
#'   \code{interactive = TRUE} or \code{"ggiraph"}, a \code{ggiraph}
#'   htmlwidget. When \code{"plotly"}, a \code{plotly} htmlwidget.
#'
#' @examples
#' \dontrun{
#' # long-format input
#' df <- tibble::tibble(
#'   barcode = rep(c("A", "B", "C"), each = 5),
#'   time    = rep(0:4, times = 3),
#'   counts  = c(100, 80, 60, 40, 20,
#'                20, 40, 60, 80, 100,
#'                 0,  0, 10, 50, 40)
#' )
#' barbac_ts_area(df, min_total_count = 0)
#'
#' # Bartender wide-format input auto-detects
#' bt <- tibble::tibble(
#'   Cluster.ID    = c("bc1", "bc2"),
#'   Center        = c("ACGT", "ACGA"),
#'   Cluster.Score = c(0.9, 0.8),
#'   time_point_1  = c(100, 20),
#'   time_point_2  = c(80,  40),
#'   time_point_3  = c(60,  60)
#' )
#' barbac_ts_area(bt, min_total_count = 0)
#' }
#'
#' @export
#' @importFrom dplyr filter mutate arrange left_join group_by ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom readr read_csv
#' @importFrom PNWColors pnw_palette
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes geom_area scale_fill_manual guides labs
#'   scale_x_continuous
barbac_ts_area <- function(data,
                           id_col           = "barcode",
                           time_col         = "time",
                           count_col        = "counts",
                           wide_time_prefix = "time_point_",
                           min_total_count  = 10,
                           include_late     = TRUE,
                           fill_missing     = c("epsilon", "zero"),
                           epsilon          = 1e-6,
                           palette          = NULL,
                           time_zero_shift  = FALSE,
                           x_breaks         = NULL,
                           title            = NULL,
                           x_lab            = "Time",
                           y_lab            = "Barcode frequency",
                           show_legend      = FALSE,
                           theme            = theme_barbac(),
                           interactive      = FALSE) {

  fill_missing <- match.arg(fill_missing)

  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) stop("File not found: ", data)
    data <- readr::read_csv(data, show_col_types = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame/tibble or a path to a CSV.")
  }

  cols <- names(data)
  is_wide <- "Cluster.ID" %in% cols &&
             any(startsWith(cols, wide_time_prefix))
  is_long <- all(c(id_col, time_col, count_col) %in% cols)

  if (!is_wide && !is_long) {
    stop(
      "Unrecognised input layout. Expected one of:\n",
      "  * Long format with columns '", id_col, "', '", time_col,
      "', '", count_col, "'.\n",
      "  * Bartender-wide with column 'Cluster.ID' plus one column per\n",
      "    timepoint starting with '", wide_time_prefix,
      "' (e.g. '", wide_time_prefix, "1', '", wide_time_prefix, "2', ...).\n",
      "Got columns: ", paste(cols, collapse = ", ")
    )
  }

  if (is_wide) {
    drop_cols <- intersect(c("Center", "Cluster.Score"), cols)
    if (length(drop_cols)) {
      data <- data[, setdiff(names(data), drop_cols), drop = FALSE]
    }
    names(data)[names(data) == "Cluster.ID"] <- id_col
    time_cols <- grep(paste0("^", wide_time_prefix), names(data), value = TRUE)
    data <- tidyr::pivot_longer(
      data,
      cols      = dplyr::all_of(time_cols),
      names_to  = time_col,
      values_to = count_col
    )
    data[[time_col]] <- as.numeric(
      sub(paste0("^", wide_time_prefix), "", data[[time_col]])
    )
  }

  data[[time_col]]  <- as.numeric(data[[time_col]])
  data[[count_col]] <- suppressWarnings(as.numeric(data[[count_col]]))
  data <- data[!is.na(data[[time_col]]), , drop = FALSE]

  if (time_zero_shift) {
    data[[time_col]][data[[time_col]] == 1] <- 0
  }

  # Guarantee one row per (barcode, time). If the shift or the input
  # itself contains duplicate (barcode, time) pairs, sum their counts;
  # otherwise geom_area with group = barcode draws crossing polygons.
  data <- stats::aggregate(
    data[[count_col]],
    by = list(data[[id_col]], data[[time_col]]),
    FUN = sum, na.rm = TRUE
  )
  names(data) <- c(id_col, time_col, count_col)

  fill_val <- if (fill_missing == "zero") 0 else epsilon

  if (include_late) {
    grid <- expand.grid(
      unique(data[[id_col]]),
      unique(data[[time_col]]),
      KEEP.OUT.ATTRS   = FALSE,
      stringsAsFactors = FALSE
    )
    names(grid) <- c(id_col, time_col)
    data <- dplyr::left_join(grid, data, by = c(id_col, time_col))
    data[[count_col]][is.na(data[[count_col]])] <- fill_val
  } else {
    min_time <- min(data[[time_col]], na.rm = TRUE)
    seen_at_start <- unique(data[[id_col]][data[[time_col]] == min_time])
    data <- data[data[[id_col]] %in% seen_at_start, , drop = FALSE]
    data[[count_col]][is.na(data[[count_col]])] <- fill_val
  }

  tot <- tapply(data[[count_col]], data[[id_col]], sum, na.rm = TRUE)
  keep_ids <- names(tot)[tot >= min_total_count]
  data <- data[data[[id_col]] %in% keep_ids, , drop = FALSE]
  if (!nrow(data)) {
    stop("No lineages left after filtering with min_total_count = ",
         min_total_count, ".")
  }

  time_tot <- tapply(data[[count_col]], data[[time_col]], sum, na.rm = TRUE)
  freq <- data[[count_col]] / time_tot[as.character(data[[time_col]])]
  freq[is.nan(freq) | is.na(freq)] <- 0
  data[[".freq"]] <- freq

  data[[id_col]] <- factor(data[[id_col]])
  n_ids <- nlevels(data[[id_col]])

  if (is.null(palette)) {
    base <- PNWColors::pnw_palette("Sailboat", n = n_ids, type = "continuous")
    colours <- grDevices::colorRampPalette(base)(n_ids)
  } else if (length(palette) < n_ids) {
    colours <- grDevices::colorRampPalette(palette)(n_ids)
  } else {
    colours <- palette[seq_len(n_ids)]
  }
  names(colours) <- levels(data[[id_col]])

  if (is.null(x_breaks)) {
    x_range <- range(data[[time_col]], na.rm = TRUE)
    x_breaks <- seq(x_range[1L], x_range[2L], by = 1)
  }

  data <- data[order(data[[time_col]]), , drop = FALSE]

  # Normalise interactive: FALSE | TRUE | "ggiraph" | "plotly".
  backend <- .resolve_interactive(interactive)

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x     = .data[[time_col]],
                 y     = .data[[".freq"]],
                 group = .data[[id_col]],
                 fill  = .data[[id_col]])
  )

  if (backend == "ggiraph") {
    # Precompute a per-row tooltip; ggiraph maps it via aes(tooltip=)
    # and renders each polygon segment with its own hover text.
    data[[".tooltip"]] <- sprintf(
      "%s\ntime %s: %.3f",
      as.character(data[[id_col]]),
      format(data[[time_col]]),
      data[[".freq"]]
    )
    p$data <- data                             # replace with tooltip-enriched frame
    p <- p +
      ggiraph::geom_area_interactive(
        ggplot2::aes(tooltip = .data[[".tooltip"]],
                     data_id = .data[[id_col]]),
        position = "stack", colour = NA
      )
  } else {
    # Static or plotly path: standard geom_area.
    # colour = NA removes the per-polygon outline so high-diversity
    # plots (hundreds+ lineages) render as smooth bands rather than
    # a mesh of dark borders.
    p <- p + ggplot2::geom_area(position = "stack", colour = NA)
  }

  p <- p +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = title, x = x_lab, y = y_lab, fill = id_col)

  if (!show_legend) p <- p + ggplot2::guides(fill = "none")
  if (!is.null(theme)) p <- p + theme

  switch(backend,
    "static"  = p,
    "ggiraph" = ggiraph::girafe(
      ggobj    = p,
      width_svg  = 9,
      height_svg = 5.5,
      options    = list(
        ggiraph::opts_hover(css = "stroke:#333333;stroke-width:1px;"),
        ggiraph::opts_hover_inv(css = "opacity:0.35;")
      )
    ),
    "plotly"  = plotly::ggplotly(p, tooltip = c("fill", "x", "y"))
  )
}


# Internal: normalise the user-facing `interactive` argument.
.resolve_interactive <- function(x) {
  if (isFALSE(x))                  return("static")
  if (isTRUE(x))                   x <- "ggiraph"
  if (!is.character(x) || length(x) != 1L)
    stop("`interactive` must be FALSE, TRUE, 'ggiraph', or 'plotly'.",
         call. = FALSE)
  x <- tolower(x)
  if (!x %in% c("ggiraph", "plotly"))
    stop("`interactive` must be FALSE, TRUE, 'ggiraph', or 'plotly'.",
         call. = FALSE)
  pkg <- x
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("interactive = '%s' requires the '%s' package; ", x, pkg),
      sprintf("install it with install.packages(\"%s\").", pkg),
      call. = FALSE
    )
  }
  x
}
