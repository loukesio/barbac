#' @title Fast Centroid-Based Sequence Clustering
#'
#' @description Abundance-ranked centroid clustering with likelihood-based
#'   best-parent selection and a distance-aware count-ratio merge guard.
#'   Supports both Levenshtein (handles indels) and Hamming distances.
#'
#' @param input_path Character string or data.frame.
#' @param distance Numeric. Maximum edit distance. Default: 3.
#' @param method Character string. "lv" (Levenshtein, default) or "hamming".
#' @param barcode_col Character string. Default: "barcode".
#' @param counts_col Character string. Default: "counts".
#' @param output_dir Character string or NULL. Default: NULL.
#' @param file_pattern Character string. Default: "\\.csv$".
#' @param verbose Logical. Default: TRUE.
#' @param use_cpp Logical. Default: TRUE.
#' @param use_kmer_filter Logical. Default: TRUE.
#' @param kmer_size Integer. Seed size for LV index. Default: 5.
#' @param min_shared_kmers Integer. Kept for API compatibility. Default: 2.
#' @param merge_ratio Numeric. Base count-ratio for the distance-aware merge
#'   guard. Effective ratio increases with distance. Default: 20.
#' @param error_rate Numeric. Approximate per-base error rate for likelihood
#'   scoring. Default: 0.005.
#' @param q Integer. Q-gram size for qgram/jaccard/cosine. Default: 2.
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   cluster_id, central_barcode, all_barcodes, all_counts, sum_counts.
#'
#' @export
#' @importFrom dplyr arrange desc
#' @importFrom tibble tibble
#' @importFrom readr read_csv write_csv
#' @importFrom stringdist stringdist
#' @importFrom rlang sym
#' @useDynLib barbac, .registration = TRUE
super_cluster2 <- function(input_path,
                           distance         = 3,
                           method           = c("lv", "hamming", "osa", "dl", "lcs",
                                                "qgram", "cosine", "jaccard", "jw", "soundex"),
                           barcode_col      = "barcode",
                           counts_col       = "counts",
                           output_dir       = NULL,
                           file_pattern     = "\\.csv$",
                           verbose          = TRUE,
                           use_cpp          = TRUE,
                           use_kmer_filter  = TRUE,
                           kmer_size        = 5L,
                           min_shared_kmers = 2L,
                           merge_ratio      = 20.0,
                           error_rate       = 0.005,
                           q                = 2) {
  
  method        <- match.arg(method)
  use_cpp_final <- use_cpp && (method %in% c("lv", "hamming"))
  
  if (is.data.frame(input_path)) {
    return(.process_df(input_path, distance, method, barcode_col, counts_col,
                       output_dir, verbose, use_cpp_final, use_kmer_filter,
                       kmer_size, min_shared_kmers, merge_ratio, error_rate, q))
  }
  
  if (!is.character(input_path))
    stop("input_path must be a file path or data.frame")
  if (!file.exists(input_path))
    stop("Input path does not exist: ", input_path)
  
  if (file.info(input_path)$isdir) {
    .process_dir(input_path, distance, method, barcode_col, counts_col,
                 output_dir, file_pattern, verbose, use_cpp_final,
                 use_kmer_filter, kmer_size, min_shared_kmers,
                 merge_ratio, error_rate, q)
  } else {
    .process_file(input_path, distance, method, barcode_col, counts_col,
                  output_dir, verbose, use_cpp_final, use_kmer_filter,
                  kmer_size, min_shared_kmers, merge_ratio, error_rate, q)
  }
}

# =============================================================================
# Internal: process a data.frame
# =============================================================================
#' @keywords internal
#' @noRd
.process_df <- function(data, distance, method, barcode_col, counts_col,
                        output_dir, verbose, use_cpp_final, use_kmer_filter,
                        kmer_size, min_shared_kmers, merge_ratio, error_rate, q) {
  
  if (!all(c(barcode_col, counts_col) %in% colnames(data)))
    stop(sprintf("Columns '%s' and/or '%s' not found. Available: %s",
                 barcode_col, counts_col,
                 paste(colnames(data), collapse = ", ")))
  
  n_before <- nrow(data)
  data     <- data[!is.na(data[[barcode_col]]) & !is.na(data[[counts_col]]), ]
  
  data[[barcode_col]] <- as.character(data[[barcode_col]])
  data[[counts_col]]  <- as.integer(data[[counts_col]])
  
  data <- dplyr::arrange(data, dplyr::desc(!!rlang::sym(counts_col)))
  
  mean_len <- mean(nchar(data[[barcode_col]]))
  
  if (verbose) {
    message("========================================")
    message("barbac: Abundance-ranked centroid clustering")
    message("  Implementation   : ",
            ifelse(use_cpp_final, "C++ (optimized)", "R (stringdist)"))
    message("  K-mer filter     : ",
            ifelse(use_kmer_filter && use_cpp_final, "ON", "OFF"))
    message("  Method           : ", method)
    message("  Distance         : ", distance)
    message("  kmer_size        : ", kmer_size)
    message("  min_shared_kmers : ", min_shared_kmers)
    message("  merge_ratio      : ", merge_ratio, " (distance-aware)")
    message("  error_rate       : ", error_rate)
    if (n_before > nrow(data))
      message("  Removed NAs      : ", n_before - nrow(data), " rows")
    message("  Sequences        : ", format(nrow(data), big.mark = ","))
    message("  Mean length      : ", round(mean_len, 1), " bp")
    message("  Top sequence     : ", data[[barcode_col]][1],
            " (count: ", format(data[[counts_col]][1], big.mark = ","), ")")
    message("Running...")
  }
  
  t0 <- Sys.time()
  
  if (use_cpp_final) {
    cpp <- barbac_cpp_centroid_cluster_optimized(
      barcodes         = data[[barcode_col]],
      counts           = data[[counts_col]],
      max_distance     = distance,
      method           = method,
      kmer_size        = kmer_size,
      min_shared_kmers = min_shared_kmers,
      use_kmer_filter  = use_kmer_filter,
      merge_ratio      = merge_ratio,
      error_rate       = error_rate,
      verbose          = verbose
    )
    result <- tibble::tibble(
      cluster_id      = cpp$cluster_id,
      central_barcode = cpp$central_barcode,
      all_barcodes    = cpp$all_barcodes,
      all_counts      = cpp$all_counts,
      sum_counts      = cpp$sum_counts
    )
    attr(result, "blocked_by_dist") <- cpp$blocked_by_dist
    attr(result, "build_id")        <- cpp$build_id
    attr(result, "method")          <- cpp$method
  } else {
    result <- .cluster_stringdist(data[[barcode_col]], data[[counts_col]],
                                  distance, method, q)
  }
  
  t1    <- Sys.time()
  tsecs <- as.numeric(difftime(t1, t0, units = "secs"))
  
  if (verbose) {
    message("\u2713 Done in ", round(tsecs, 1), "s (",
            round(tsecs / 60, 1), " min)")
    message("  Clusters found   : ", format(nrow(result), big.mark = ","))
    message("  Compression      : ",
            round(nrow(data) / nrow(result), 2), ":1")
    message("  Rate             : ",
            format(round(nrow(data) / tsecs), big.mark = ","), " seq/s")
    if (use_cpp_final) {
      message("  Build ID         : ", attr(result, "build_id"))
      bd <- attr(result, "blocked_by_dist")
      if (!is.null(bd)) {
        for (d in seq_along(bd)) {
          if (bd[d] > 0) {
            message("  Blocked at d=", d - 1, "  : ",
                    format(bd[d], big.mark = ","))
          }
        }
      }
    }
    message("========================================\n")
  }
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out <- result
    out$all_barcodes <- sapply(result$all_barcodes, paste, collapse = ",")
    out$all_counts   <- sapply(result$all_counts,   paste, collapse = ",")
    readr::write_csv(out, file.path(output_dir, "clustered.csv"))
    if (verbose) message("Saved to: ", output_dir)
  }
  
  result
}

# =============================================================================
# Internal: process a single CSV file
# =============================================================================
#' @keywords internal
#' @noRd
.process_file <- function(file_path, distance, method, barcode_col, counts_col,
                          output_dir, verbose, use_cpp_final, use_kmer_filter,
                          kmer_size, min_shared_kmers, merge_ratio, error_rate, q) {
  
  if (verbose) message("Reading: ", basename(file_path))
  data <- readr::read_csv(file_path, show_col_types = FALSE)
  
  if (!all(c(barcode_col, counts_col) %in% colnames(data)))
    stop("Required columns not found in: ", file_path)
  
  result <- .process_df(data, distance, method, barcode_col, counts_col,
                        NULL, verbose, use_cpp_final, use_kmer_filter,
                        kmer_size, min_shared_kmers, merge_ratio, error_rate, q)
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out_path <- file.path(output_dir,
                          paste0(tools::file_path_sans_ext(basename(file_path)),
                                 "_clustered.csv"))
    out <- result
    out$all_barcodes <- sapply(result$all_barcodes, paste, collapse = ",")
    out$all_counts   <- sapply(result$all_counts,   paste, collapse = ",")
    readr::write_csv(out, out_path)
    if (verbose) message("Saved: ", out_path)
  }
  result
}

# =============================================================================
# Internal: process a directory of CSV files
# =============================================================================
#' @keywords internal
#' @noRd
.process_dir <- function(dir_path, distance, method, barcode_col, counts_col,
                         output_dir, file_pattern, verbose, use_cpp_final,
                         use_kmer_filter, kmer_size, min_shared_kmers,
                         merge_ratio, error_rate, q) {
  
  files <- list.files(dir_path, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0)
    stop("No CSV files found matching pattern: ", file_pattern)
  if (verbose) message("Found ", length(files), " CSV files")
  
  results <- list()
  for (f in files) {
    nm <- tools::file_path_sans_ext(basename(f))
    tryCatch(
      results[[nm]] <- .process_file(f, distance, method, barcode_col,
                                     counts_col, output_dir, verbose,
                                     use_cpp_final, use_kmer_filter,
                                     kmer_size, min_shared_kmers,
                                     merge_ratio, error_rate, q),
      error = function(e) warning("Failed: ", basename(f), ": ", e$message)
    )
  }
  results
}

# =============================================================================
# R fallback (no C++)
# =============================================================================
#' @keywords internal
#' @noRd
.cluster_stringdist <- function(barcodes, counts, max_distance, method, q) {
  n <- length(barcodes)
  if (n == 0) {
    return(tibble::tibble(
      cluster_id      = character(),
      central_barcode = character(),
      all_barcodes    = list(),
      all_counts      = list(),
      sum_counts      = integer()
    ))
  }
  
  clusters <- list()
  for (i in seq_len(n)) {
    bc  <- barcodes[i]
    cnt <- counts[i]
    hit <- NULL
    for (j in seq_along(clusters)) {
      d <- if (method %in% c("qgram", "cosine", "jaccard")) {
        stringdist::stringdist(bc, clusters[[j]]$central_barcode,
                               method = method, q = q)
      } else {
        stringdist::stringdist(bc, clusters[[j]]$central_barcode,
                               method = method)
      }
      if (method == "jw") d <- 1 - d
      if (d <= max_distance) { hit <- j; break }
    }
    if (!is.null(hit)) {
      clusters[[hit]]$all_barcodes <- c(clusters[[hit]]$all_barcodes, bc)
      clusters[[hit]]$all_counts   <- c(clusters[[hit]]$all_counts, cnt)
      clusters[[hit]]$sum_counts   <- clusters[[hit]]$sum_counts + cnt
    } else {
      clusters <- c(clusters, list(list(
        cluster_id      = paste0("group", length(clusters) + 1),
        central_barcode = bc,
        all_barcodes    = bc,
        all_counts      = cnt,
        sum_counts      = cnt
      )))
    }
  }
  
  tibble::tibble(
    cluster_id      = vapply(clusters, `[[`, character(1), "cluster_id"),
    central_barcode = vapply(clusters, `[[`, character(1), "central_barcode"),
    all_barcodes    = lapply(clusters, `[[`, "all_barcodes"),
    all_counts      = lapply(clusters, `[[`, "all_counts"),
    sum_counts      = vapply(clusters, `[[`, integer(1),   "sum_counts")
  )
}