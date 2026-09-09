# Record the exact Quarto output and the source/input hashes behind it.
finalize_r_report <- function(report, root, kit, here) {
  files <- c(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE),
    file.path(here, c("build_report.R", "report_helpers.R", "finalize_report.R", "report.qmd")))
  original <- jsonlite::read_json(file.path(kit, "results", "validation.json"))$artifact_sha256
  assets <- list.files(here, pattern = "\\.(png|csv)$", recursive = TRUE)
  assets <- assets[!grepl("^qa/", assets)]
  for (name in names(original)) stopifnot(report_sha(file.path(kit, "results", name)) == original[[name]])
  receipt <- list(status = "numerical_checks_passed", language = "R", report = "report.html",
    renderer = paste("Quarto", system2(Sys.which("quarto"), "--version", stdout = TRUE)),
    checks = report$checks, complete_samples = 8, sample_method_comparisons = 48,
    package_functions = c("summarise_bam_stats", "plot_bam_stats", "barbac_xtr.stats", "cluster_stats", "barbac_ts_area", "theme_barbac"),
    upstream_results_unchanged = TRUE, browser_validation = "pending",
    historical_analysis_commit = "bf171db",
    current_base_commit = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
    source_sha256 = setNames(lapply(files, report_sha), substring(files, nchar(root) + 2)),
    input_sha256 = report$input_hashes, original_artifact_sha256 = original,
    output_sha256 = report_sha(file.path(here, "report.html")),
    figure_table_sha256 = setNames(lapply(file.path(here, assets), report_sha), assets),
    limitations = c("Publication agreement is not biological ground-truth accuracy.",
      "Pooled clustering uses all four generations; no fitness model is fitted.",
      "Shared extraction and mapping costs are excluded from clustering timings.",
      "Length and entropy histograms count distinct sequences; package labels Reads refer to UMI molecules here.",
      "The extraction summary length-bin calculation was corrected; clustering and extracted counts were not changed."))
  jsonlite::write_json(receipt, file.path(here, "validation.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
  invisible(receipt)
}
