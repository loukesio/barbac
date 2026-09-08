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
#' @param kmer_size Integer. Kept for API compatibility. Search partitions are
#'   now chosen from observed sequence information. Default: 5.
#' @param min_shared_kmers Integer. Kept for API compatibility. Default: 2.
#' @param merge_ratio Numeric. Base count-ratio for the distance-aware merge
#'   guard. Effective ratio increases with distance. Default: 20.
#' @param error_rate Numeric. Approximate per-base error rate for likelihood
#'   scoring. Default: 0.005.
#' @param tie_break Character string. How to order barcodes that share a count:
#'   \code{"sequence"} (default) orders them by the barcode itself;
#'   \code{"hash"} orders them by a salted hash using \code{tie_seed};
#'   \code{"support"} first orders by the summed counts of one-edit neighbours
#'   that are no more abundant than the barcode, then by sequence. Support
#'   uses the selected distance method and observed reads only. It is an
#'   optional evidence-based tie rule, not a guarantee of improved accuracy.
#'   All three options are deterministic and independent of input row order.
#'   Re-running the hash option across seeds measures sensitivity to arbitrary
#'   equal-count ordering; support retains sequence order when evidence ties.
#' @param tie_seed Integer. Salt for \code{tie_break = "hash"}. Default: 0.
#' @param use_design Logical. Exploit the barcode design. A library that fixes
#'   some positions and randomises others carries identity only at the random
#'   ones, so a read differing from a centroid solely at a fixed position cannot
#'   be a different barcode and is absorbed without consulting the abundance
#'   guard. The fixed positions are read off the data (a position where one base
#'   covers at least 90% of reads), so no template has to be supplied. Has no
#'   effect on fully random libraries, where every position varies.
#'   Default: FALSE.
#'
#' @param indel_model Character string. Experimental LV merge-guard exception:
#'   \code{"none"} (default) or \code{"poisson"}. The latter allows a
#'   repeated-base single insertion/deletion to pass a blocked ratio guard if
#'   its count is consistent with a Poisson error expectation (upper tail at
#'   least 0.01). Uses \code{error_rate} as a deletion-rate upper-bound proxy,
#'   and one fourth of that rate for a specific inserted base, multiplied by
#'   equivalent gap positions. This is not an estimated platform-specific
#'   error rate or posterior probability. It may merge genuine length variants;
#'   validate on independent controls. Requires the C++ LV method.
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
                           method           = c("lv", "hamming"),
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
                           tie_break        = c("sequence", "hash", "support"),
                           tie_seed         = 0L,
                           use_design       = FALSE,
                           indel_model      = c("none", "poisson")) {

  method        <- match.arg(method)
  tie_break     <- match.arg(tie_break)
  indel_model   <- match.arg(indel_model)
  if (indel_model != "none" && (!use_cpp || method != "lv"))
    stop("indel_model requires method = 'lv' and use_cpp = TRUE")
  use_cpp_final <- use_cpp && (method %in% c("lv", "hamming"))

  if (is.data.frame(input_path)) {
    return(.process_df(input_path, distance, method, barcode_col, counts_col,
                       output_dir, verbose, use_cpp_final, use_kmer_filter,
                       kmer_size, min_shared_kmers, merge_ratio, error_rate,
                       tie_break, tie_seed, use_design, indel_model))
  }

  if (!is.character(input_path))
    stop("input_path must be a file path or data.frame")
  if (!file.exists(input_path))
    stop("Input path does not exist: ", input_path)

  if (file.info(input_path)$isdir) {
    .process_dir(input_path, distance, method, barcode_col, counts_col,
                 output_dir, file_pattern, verbose, use_cpp_final,
                 use_kmer_filter, kmer_size, min_shared_kmers,
                 merge_ratio, error_rate, tie_break, tie_seed, use_design, indel_model)
  } else {
    .process_file(input_path, distance, method, barcode_col, counts_col,
                  output_dir, verbose, use_cpp_final, use_kmer_filter,
                  kmer_size, min_shared_kmers, merge_ratio, error_rate,
                  tie_break, tie_seed, use_design, indel_model)
  }
}

# =============================================================================
# Internal: process a data.frame
# =============================================================================
#' @keywords internal
#' @noRd
.process_df <- function(data, distance, method, barcode_col, counts_col,
                        output_dir, verbose, use_cpp_final, use_kmer_filter,
                        kmer_size, min_shared_kmers, merge_ratio, error_rate,
                        tie_break = "sequence", tie_seed = 0L,
                        use_design = FALSE, indel_model = "none") {

  if (!all(c(barcode_col, counts_col) %in% colnames(data)))
    stop(sprintf("Columns '%s' and/or '%s' not found. Available: %s",
                 barcode_col, counts_col,
                 paste(colnames(data), collapse = ", ")))

  n_before <- nrow(data)
  data     <- data[!is.na(data[[barcode_col]]) & !is.na(data[[counts_col]]), ]

  data[[barcode_col]] <- as.character(data[[barcode_col]])
  data[[counts_col]]  <- as.integer(data[[counts_col]])

  # Collapse exact-duplicate barcodes by summing their counts. The clustering
  # kernel treats a distance of 0 as "already the same sequence" and never
  # merges identical barcodes, so duplicate rows would otherwise be double-
  # counted as separate singleton clusters. Only rewrite the table when
  # duplicates actually exist, and preserve first-occurrence row order: the
  # abundance-ranked greedy pass is ordered explicitly below.
  n_after_na <- nrow(data)
  if (anyDuplicated(data[[barcode_col]])) {
    summed <- rowsum(data[[counts_col]], group = data[[barcode_col]],
                     reorder = FALSE)
    data <- data[!duplicated(data[[barcode_col]]), , drop = FALSE]
    data[[counts_col]] <- as.integer(summed[match(data[[barcode_col]], rownames(summed)), 1L])
  }
  n_collapsed <- n_after_na - nrow(data)

  # Hamming mode compares sequences via 2-bit DNA packing, which only supports
  # A/C/G/T barcodes of length <= 32. Such sequences would otherwise be left as
  # their own singleton clusters with no diagnostic; warn instead of silently.
  if (method == "hamming") {
    bad <- grepl("[^ACGTacgt]", data[[barcode_col]]) |
           nchar(data[[barcode_col]]) > 32L
    if (any(bad))
      warning(sum(bad), " barcode(s) contain non-ACGT characters or exceed ",
              "32 bp; Hamming mode cannot compare these and will leave each as ",
              "its own cluster. Use method = 'lv' for indel/N-containing data.",
              call. = FALSE)
  }

  # Standardise the row order so the abundance-ranked greedy pass is a pure
  # function of the input's content, not the order it happened to arrive in.
  # Breaking count ties deterministically (barcodes are unique after the dedup
  # above) makes the clustering reproducible under any row permutation from
  # upstream joins, summaries, or file merges.
  #
  # Which tied barcode is visited first is nonetheless arbitrary, and it decides
  # which of them is allowed to seed a cluster. tie_break = "hash" re-draws that
  # arbitrary order from tie_seed without reference to the bases, so repeating a
  # run across seeds measures how much of a result rests on the choice rather
  # than on the data. Every seed is itself fully reproducible.
  data <- if (tie_break == "hash") {
    dplyr::arrange(data, dplyr::desc(!!rlang::sym(counts_col)),
                   barbac_seq_order_key(!!rlang::sym(barcode_col),
                                        as.integer(tie_seed)),
                   !!rlang::sym(barcode_col))
  } else if (tie_break == "support") {
    support <- barbac_support_order_key(data[[barcode_col]],
                                        data[[counts_col]], method)
    data[order(-data[[counts_col]], -support, data[[barcode_col]],
               method = "radix"), , drop = FALSE]
  } else {
    dplyr::arrange(data, dplyr::desc(!!rlang::sym(counts_col)),
                   !!rlang::sym(barcode_col))
  }

  barcode_lengths <- nchar(data[[barcode_col]])
  mean_len <- mean(barcode_lengths)
  fixed_length <- length(unique(barcode_lengths)) == 1L
  
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
    if (n_before > n_after_na)
      message("  Removed NAs      : ", n_before - n_after_na, " rows")
    if (n_collapsed > 0)
      message("  Collapsed dups   : ", n_collapsed, " rows")
    message("  Sequences        : ", format(nrow(data), big.mark = ","))
    message("  Mean length      : ", round(mean_len, 1), " bp")
    if (method == "lv" && fixed_length) {
      message("  Speed note       : all barcodes have the same length. If the ",
              "data are known to exclude indels and shifted alignments, ",
              "method = \"hamming\" is usually much faster.")
    }
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
      verbose          = verbose,
      use_design       = use_design,
      use_indel_model  = indel_model == "poisson"
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
                                  distance, method)
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
                          kmer_size, min_shared_kmers, merge_ratio, error_rate,
                          tie_break = "sequence", tie_seed = 0L,
                          use_design = FALSE, indel_model = "none") {

  if (verbose) message("Reading: ", basename(file_path))
  data <- readr::read_csv(file_path, show_col_types = FALSE)

  if (!all(c(barcode_col, counts_col) %in% colnames(data)))
    stop("Required columns not found in: ", file_path)

  result <- .process_df(data, distance, method, barcode_col, counts_col,
                        NULL, verbose, use_cpp_final, use_kmer_filter,
                        kmer_size, min_shared_kmers, merge_ratio, error_rate,
                        tie_break, tie_seed, use_design, indel_model)
  
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
                         merge_ratio, error_rate,
                         tie_break = "sequence", tie_seed = 0L,
                         use_design = FALSE, indel_model = "none") {
  
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
                                     merge_ratio, error_rate,
                                     tie_break, tie_seed, use_design, indel_model),
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
.cluster_stringdist <- function(barcodes, counts, max_distance, method) {
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
      d <- stringdist::stringdist(bc, clusters[[j]]$central_barcode,
                                  method = method)
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
