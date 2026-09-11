# Report adapters only; clustering and extraction inputs are never modified.
report_labels <- c(barbac_lv = "barbac LV + Poisson", barbac_hamming = "barbac Hamming",
  shepherd = "Shepherd (e = 0.005)", starcode_sphere = "Starcode sphere",
  starcode_mp = "Starcode message passing", bartender = "Bartender")
report_blue <- "#2C4A63"
read_report_csv <- function(path) as.data.frame(readr::read_csv(path, show_col_types = FALSE,
  progress = FALSE, num_threads = 1))
report_sha <- function(path) digest::digest(file = path, algo = "sha256")

report_plotly <- function(plot, tooltip) {
  widget <- plotly::ggplotly(plot, tooltip = tooltip)
  widget$width <- "100%"
  widget$height <- 480
  widget$x$layout$width <- NULL
  widget$x$layout$autosize <- TRUE
  widget$sizingPolicy$knitr$figure <- FALSE
  widget$sizingPolicy$knitr$defaultWidth <- "100%"
  widget$sizingPolicy$knitr$defaultHeight <- 480
  widget
}

positive_cluster_stats <- function(counts) {
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)))
  barbac::cluster_stats(data.frame(sum_counts = counts[counts > 0]), verbose = FALSE)
}

composition_input <- function(series, ids, labels) {
  stopifnot(!anyDuplicated(series[c("Barcode", "generation")]), all(series$counts >= 0))
  group <- ifelse(series$Barcode %in% ids, labels[series$Barcode], "Other")
  result <- stats::aggregate(series$counts, list(barcode = unname(group), time = series$generation), sum)
  names(result)[3] <- "counts"
  grid <- expand.grid(barcode = c(unname(labels[ids]), "Other"), time = sort(unique(series$generation)),
                      stringsAsFactors = FALSE)
  result <- merge(grid, result, all.x = TRUE, sort = FALSE)
  result$counts[is.na(result$counts)] <- 0
  expected <- tapply(series$counts, series$generation, sum)
  actual <- tapply(result$counts, result$time, sum)
  stopifnot(isTRUE(all.equal(actual, expected, check.attributes = FALSE)))
  result
}

report_table <- function(data, caption = NULL, digits = 2, page_length = 8, filter = "none") {
  table <- DT::datatable(data, rownames = FALSE, caption = caption, filter = filter,
    extensions = "Buttons", class = "stripe hover compact", escape = TRUE,
    options = list(pageLength = page_length, lengthMenu = c(8, 20, 50, 100),
      scrollX = TRUE, deferRender = TRUE, dom = "Blfrtip",
      buttons = list(list(extend = "csvHtml5", text = "Download CSV", title = "barbac_report",
        exportOptions = list(orthogonal = "export", modifier = list(search = "applied", page = "all"))))))
  numeric <- names(data)[vapply(data, is.numeric, logical(1))]
  integers <- numeric[vapply(data[numeric], function(x) all(is.na(x) | x == floor(x)), logical(1))]
  if (length(integers)) table <- DT::formatRound(table, integers, digits = 0)
  if (length(setdiff(numeric, integers))) table <- DT::formatRound(table, setdiff(numeric, integers), digits = digits)
  table
}

prepare_r_report <- function(root, kit, work, output) {
  results <- file.path(kit, "results")
  validation <- jsonlite::read_json(file.path(results, "validation.json"), simplifyVector = TRUE)
  stopifnot(validation$status == "numerical_checks_passed", validation$visual_inspection == "passed")
  # These are the immutable artifacts of the previous analysis; its historical
  # source hashes are not rewritten to pretend that later QC fixes ran earlier.
  for (name in names(validation$artifact_sha256)) {
    stopifnot(report_sha(file.path(results, name)) == validation$artifact_sha256[[name]])
  }
  manifest <- read.delim(file.path(kit, "samples.tsv"), stringsAsFactors = FALSE)
  extraction <- read_report_csv(file.path(results, "extraction_summary.csv"))
  agreement <- read_report_csv(file.path(results, "sample_method_agreement.csv"))
  timings <- read_report_csv(file.path(results, "clustering_runs.csv"))
  original_summary <- read_report_csv(file.path(results, "method_summary.csv"))
  provenance <- jsonlite::read_json(file.path(results, "provenance.json"), simplifyVector = TRUE)
  dir.create(file.path(output, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output, "tables"), recursive = TRUE, showWarnings = FALSE)
  input_hashes <- list(); mapping <- list(); panels <- list(); extraction_stats <- list()
  components <- c(diverse = "BC2 · diverse barcode", environment = "BC1 · environment barcode")
  barbac::use_barbac_env()
  stopifnot(nzchar(Sys.which("samtools")))
  for (i in seq_len(nrow(manifest))) {
    sample <- manifest$sample[i]
    message("R report: mapping and extraction summaries for ", sample)
    directory <- file.path(work, "mapping", sample, "merged", "bam")
    bam <- file.path(directory, paste0(manifest$run[i], "_ANC.assembled_sorted.bam"))
    stats <- barbac::summarise_bam_stats(directory)
    stopifnot(nrow(stats) == 1L, all(is.finite(unlist(stats[c("mapped", "unmapped")]))))
    receipt <- jsonlite::read_json(file.path(work, "mapping", sample, "mapping.json"), simplifyVector = TRUE)
    if (!is.null(receipt$stats)) {
      stopifnot(stats$mapped == receipt$stats$mapped, stats$unmapped == receipt$stats$unmapped)
    }
    stats$sample <- sample
    stats$input_pairs <- extraction$total_pairs[match(sample, extraction$sample)]
    stats$merged_alignments <- stats$mapped + stats$unmapped
    stats$mapped_percent_merged <- 100 * stats$mapped / stats$merged_alignments
    mapping[[sample]] <- stats
    input_hashes[[paste0(sample, "/mapping.json")]] <- report_sha(file.path(work, "mapping", sample, "mapping.json"))
    extracted <- file.path(work, "extracted", sample)
    receipt <- jsonlite::read_json(file.path(extracted, "extraction.json"), simplifyVector = TRUE)
    stopifnot(unname(tools::md5sum(bam)) == receipt$bam_extraction$bam_md5)
    input_hashes[[paste0(sample, "/extraction.json")]] <- report_sha(file.path(extracted, "extraction.json"))
    for (component in names(components)) {
      name <- paste0(component, "_barbac.csv")
      input <- file.path(extracted, name)
      stopifnot(report_sha(input) == receipt$outputs[[name]])
      input_hashes[[paste0(sample, "/", name)]] <- report_sha(input)
      data <- read_report_csv(input)
      stopifnot(!anyDuplicated(data$barcode), all(data$counts > 0),
                all(nchar(data$barcode) == data$barcode_length),
                sum(data$counts) == receipt$barbac_input_molecules)
      panel <- barbac::barbac_xtr.stats(data, barcode_length = c(24, 28),
        fill_color = report_blue, verbose = FALSE)
      # The existing package supplies the three histograms and length-bin table.
      figure <- file.path(output, "figures", paste0(sample, "_", component, ".png"))
      ggplot2::ggsave(figure, panel, width = 11, height = 7, dpi = 125, device = ragg::agg_png)
      panels[[paste(sample, component)]] <- list(sample = sample, component = unname(components[component]),
                                                path = figure)
      extraction_stats[[paste(sample, component)]] <- data.frame(sample = sample, component = component,
        distinct_sequences = nrow(data), molecules = sum(data$counts), min_length = min(data$barcode_length),
        max_length = max(data$barcode_length), mean_length_unique = mean(data$barcode_length),
        median_abundance = median(data$counts), mean_abundance = mean(data$counts))
    }
  }
  mapping <- do.call(rbind, mapping); rownames(mapping) <- NULL
  extraction_stats <- do.call(rbind, extraction_stats); rownames(extraction_stats) <- NULL
  series <- list(); recalculated <- list(); cluster_pairs <- list(); cluster_components <- list()
  explorer <- list()
  for (method in names(report_labels)) {
    message("R report: validating ", method)
    x <- read_report_csv(file.path(results, paste0("time_series_", method, ".csv.gz")))
    stopifnot(!anyDuplicated(x[c("Barcode", "sample")]), all(x$counts >= 0))
    for (sample in manifest$sample) {
      frame <- x[x$sample == sample, ]; reference <- frame[frame$is_published_pair, ]
      stopifnot(nrow(reference) == 2314, length(unique(frame$assigned_molecules)) == 1,
                sum(frame$counts) == frame$assigned_molecules[1],
                abs(sum(frame$frequency_assigned) - 1) < 1e-10)
      expected <- agreement[agreement$method == method & agreement$sample == sample, ]
      rho <- stats::cor(reference$counts, reference$published_counts, method = "spearman")
      mass <- 100 * sum(reference$counts) / frame$input_molecules[1]
      tv <- sum(abs(reference$counts / sum(reference$counts) -
                      reference$published_counts / sum(reference$published_counts))) / 2
      stopifnot(abs(rho - expected$spearman_published_pairs) < 1e-12,
        abs(mass - expected$published_pair_mass_percent) < 1e-10,
        abs(tv - expected$total_variation_published_set) < 1e-12)
      recalculated[[paste(method, sample)]] <- data.frame(method = method, sample = sample,
        spearman = rho, published_mass_percent = mass, total_variation = tv)
      if (method %in% c("barbac_lv", "barbac_hamming")) {
        s <- positive_cluster_stats(frame$counts)
        cluster_pairs[[paste(method, sample)]] <- cbind(method = method, sample = sample, s)
      }
    }
    if (method %in% c("barbac_lv", "barbac_hamming")) {
      series[[method]] <- x
      for (component in names(components)) {
        path <- file.path(work, "clustering", "repeat_1", component, method, "centroids.csv")
        receipt <- jsonlite::read_json(file.path(dirname(path), "run.json"), simplifyVector = TRUE)
        stopifnot(report_sha(path) == receipt$output_sha256[["centroids.csv"]])
        s <- barbac::cluster_stats(path, verbose = FALSE)
        stopifnot(s$total_reads == sum(extraction$barbac_input_molecules))
        cluster_components[[paste(method, component)]] <- cbind(method = method, component = component, s)
        input_hashes[[paste(method, component, "centroids.csv", sep = "/")]] <- report_sha(path)
      }
      wide <- tidyr::pivot_wider(x[c("Barcode", "sample", "counts")], names_from = "sample", values_from = "counts")
      wide <- as.data.frame(wide)
      count_cols <- match(manifest$sample, names(wide))
      wide$pooled_molecules <- rowSums(wide[count_cols])
      wide$is_published_pair <- wide$Barcode %in% x$Barcode[x$is_published_pair]
      wide <- wide[order(-wide$pooled_molecules, wide$Barcode), c("Barcode", "is_published_pair", "pooled_molecules", manifest$sample)]
      names(wide)[-(1:3)] <- paste0(manifest$replicate, " · generation ", manifest$generation)
      explorer[[method]] <- wide
    }
  }
  recalculated <- do.call(rbind, recalculated); rownames(recalculated) <- NULL
  summary <- do.call(rbind, lapply(names(report_labels), function(method) {
    a <- recalculated[recalculated$method == method, ]
    time <- aggregate(workflow_seconds ~ `repeat`, data = timings[timings$method == method, ], sum)$workflow_seconds
    stopifnot(length(time) == 3)
    data.frame(method = method, spearman_median = median(a$spearman),
      published_mass_percent_median = median(a$published_mass_percent),
      clustering_seconds_median = median(time), clustering_seconds_min = min(time),
      clustering_seconds_max = max(time), unassigned_molecules = sum(agreement$unassigned_molecules[agreement$method == method]))
  }))
  for (column in setdiff(names(summary), "method")) {
    stopifnot(isTRUE(all.equal(summary[[column]], original_summary[[column]][match(summary$method, original_summary$method)],
                              tolerance = 1e-10, check.attributes = FALSE)))
  }
  cluster_components <- do.call(rbind, cluster_components); rownames(cluster_components) <- NULL
  cluster_pairs <- do.call(rbind, cluster_pairs); rownames(cluster_pairs) <- NULL
  lv <- series$barbac_lv
  pooled <- sort(tapply(lv$counts, lv$Barcode, sum), decreasing = TRUE)
  top <- names(pooled)[1:4]
  lookup <- read_report_csv(file.path(results, "figure_lineages.csv"))
  labels <- setNames(lookup$display_id, lookup$Barcode)
  stopifnot(all(top %in% names(labels)))
  palette <- setNames(c("#2C4A63", "#B07D22", "#7B8245", "#BD668A", "#D4D4D4"), c(unname(labels[top]), "Other"))
  # barbac_ts_area orders its factor levels alphabetically; supply matching colors.
  palette <- palette[sort(names(palette))]
  compositions <- list()
  for (method in names(series)) for (replicate in c("R1", "R2")) {
    input <- composition_input(series[[method]][series[[method]]$replicate == replicate, ], top, labels)
    args <- list(data = input, min_total_count = 0, include_late = TRUE,
      fill_missing = "zero", time_zero_shift = FALSE, x_breaks = c(8, 16, 24, 40),
      palette = unname(palette), show_legend = TRUE,
      x_lab = "Generation", y_lab = "Fraction of all assigned molecules",
      title = paste(report_labels[method], replicate, sep = " · "), theme = barbac::theme_barbac(base_size = 12, family = "Arial"))
    static <- do.call(barbac::barbac_ts_area, args)
    stopifnot(all(abs(tapply(static$data$.freq, static$data$time, sum) - 1) < 1e-12),
      all(static$data$.freq[static$data$counts == 0] == 0),
      setequal(unique(static$data$time), c(8, 16, 24, 40)))
    file <- file.path(output, "figures", paste0(method, "_", replicate, "_composition.png"))
    ggplot2::ggsave(file, static, width = 10, height = 5.5, dpi = 130, device = ragg::agg_png)
    widget <- do.call(barbac::barbac_ts_area, c(args, list(interactive = "ggiraph")))
    widget <- ggiraph::girafe_options(widget, ggiraph::opts_sizing(rescale = TRUE),
      ggiraph::opts_zoom(min = 1, max = 5), ggiraph::opts_toolbar(saveaspng = TRUE))
    compositions[[paste(method, replicate)]] <- list(widget = widget, path = file,
      method = method, replicate = replicate, input = input)
  }
  tables <- list(mapping = mapping, extraction = extraction_stats, component_clusters = cluster_components,
    paired_clusters = cluster_pairs, agreement = recalculated, methods = summary)
  for (name in names(tables)) readr::write_csv(tables[[name]], file.path(output, "tables", paste0(name, ".csv")))
  report <- list(manifest = manifest, mapping = mapping, extraction = extraction, extraction_stats = extraction_stats,
    panels = panels, cluster_components = cluster_components, cluster_pairs = cluster_pairs,
    summary = summary, agreement = recalculated, original_agreement = agreement, explorer = explorer,
    compositions = compositions, lookup = lookup, provenance = provenance,
    unmatched = jsonlite::read_json(file.path(results, "unmatched_diagnosis.json"), simplifyVector = TRUE),
    input_hashes = input_hashes, output = output,
    checks = list(source_artifact_hashes_verified = TRUE, mapping_counts_recomputed_with_barbac = TRUE,
      extraction_component_molecule_totals_match = TRUE, correlations_recomputed_in_R = TRUE,
      coverage_and_total_variation_recomputed_in_R = TRUE, timing_summary_recomputed_in_R = TRUE,
      zero_counts_excluded_from_cluster_statistics = TRUE, composition_totals_preserved = TRUE,
      four_actual_generations_preserved = TRUE, no_pseudocounts = TRUE))
  report
}
