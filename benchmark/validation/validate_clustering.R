# =============================================================================
# 02_validation_functions.R
# Clustering validation against ground truth
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringdist)
library(RColorBrewer)

# -----------------------------------------------------------------------------
# Core validation function
# -----------------------------------------------------------------------------

#' Evaluate clustering results against known ground truth
#' 
#' @param cluster_result Output from super_cluster2()
#' @param barcode_data Original data with parent_barcode column (ground truth)
#' @param parent_barcodes Vector of true parent sequences
#' 
#' @return List with purity, completeness, parent_recovery_rate, f1_score
#' @export
evaluate_single_clustering <- function(cluster_result, barcode_data, parent_barcodes) {
  
  barcode_to_parent <- barcode_data %>%
    dplyr::select(barcode, parent_barcode) %>%
    dplyr::distinct() %>%
    tibble::deframe()
  
  parent_to_id <- setNames(seq_along(parent_barcodes), parent_barcodes)
  
  # PURITY
  cluster_purities <- numeric(nrow(cluster_result))
  cluster_dominant_parent <- character(nrow(cluster_result))
  
  for (i in 1:nrow(cluster_result)) {
    cluster_barcodes <- cluster_result$all_barcodes[[i]]
    cluster_counts <- cluster_result$all_counts[[i]]
    
    parents_in_cluster <- sapply(cluster_barcodes, function(bc) {
      barcode_to_parent[bc]
    })
    
    parent_counts <- tapply(cluster_counts, parents_in_cluster, sum)
    total_counts <- sum(parent_counts)
    max_parent_count <- max(parent_counts)
    cluster_purities[i] <- max_parent_count / total_counts
    cluster_dominant_parent[i] <- names(which.max(parent_counts))
  }
  
  cluster_sizes <- sapply(cluster_result$all_counts, sum)
  overall_purity <- sum(cluster_purities * cluster_sizes) / sum(cluster_sizes)
  
  # COMPLETENESS
  parent_completeness <- numeric(length(parent_barcodes))
  
  for (p in seq_along(parent_barcodes)) {
    parent_seq <- parent_barcodes[p]
    
    offspring_data <- barcode_data %>%
      dplyr::filter(parent_barcode == parent_seq)
    
    if (nrow(offspring_data) == 0) {
      parent_completeness[p] <- NA
      next
    }
    
    offspring_barcodes <- offspring_data$barcode
    offspring_counts <- offspring_data$counts
    
    cluster_assignments <- numeric(length(offspring_barcodes))
    
    for (j in seq_along(offspring_barcodes)) {
      bc <- offspring_barcodes[j]
      for (k in 1:nrow(cluster_result)) {
        if (bc %in% cluster_result$all_barcodes[[k]]) {
          cluster_assignments[j] <- k
          break
        }
      }
    }
    
    cluster_counts_per_parent <- tapply(offspring_counts, cluster_assignments, sum)
    total_offspring_counts <- sum(offspring_counts)
    max_cluster_count <- max(cluster_counts_per_parent)
    
    parent_completeness[p] <- max_cluster_count / total_offspring_counts
  }
  
  overall_completeness <- mean(parent_completeness, na.rm = TRUE)
  
  # PARENT RECOVERY
  centroids <- cluster_result$central_barcode
  parents_recovered <- 0
  recovery_details <- list()
  
  for (p in seq_along(parent_barcodes)) {
    parent_seq <- parent_barcodes[p]
    
    distances_to_centroids <- stringdist::stringdist(parent_seq, centroids, method = "lv")
    min_distance <- min(distances_to_centroids)
    closest_centroid <- centroids[which.min(distances_to_centroids)]
    
    recovered <- min_distance <= 2
    if (recovered) parents_recovered <- parents_recovered + 1
    
    recovery_details[[p]] <- list(
      parent = parent_seq,
      closest_centroid = closest_centroid,
      distance = min_distance,
      recovered = recovered
    )
  }
  
  parent_recovery_rate <- parents_recovered / length(parent_barcodes)
  
  # F1 SCORE
  f1_score <- 2 * (overall_purity * overall_completeness) / 
    (overall_purity + overall_completeness)
  
  return(list(
    purity = overall_purity,
    completeness = overall_completeness,
    parent_recovery_rate = parent_recovery_rate,
    f1_score = f1_score,
    cluster_purities = cluster_purities,
    cluster_dominant_parent = cluster_dominant_parent,
    parent_completeness = parent_completeness,
    recovery_details = recovery_details
  ))
}

# -----------------------------------------------------------------------------
# Cluster and validate all time points
# -----------------------------------------------------------------------------

#' Cluster and validate all time points from a simulation
#' 
#' @param sim_result Output from simulate_timeseries_*()
#' @param distance Distance threshold for super_cluster2()
#' @param method Distance method (default "lv")
#' @param verbose Print progress
#' 
#' @return List with clustering results and validation metrics per time point
#' @export
cluster_and_validate_timeseries <- function(sim_result, 
                                            distance = 2,
                                            method = "lv",
                                            verbose = TRUE) {
  
  parent_barcodes <- sim_result$parent_barcodes
  timepoints <- names(sim_result$timepoint_data)
  
  if (verbose) {
    cat("\n")
    cat("================================================================\n")
    cat("     CLUSTERING AND VALIDATION ACROSS TIME SERIES              \n")
    cat("================================================================\n")
    cat(sprintf("  Model: %s\n", sim_result$model))
    cat(sprintf("  Distance threshold: %d\n", distance))
    cat(sprintf("  Time points: %d\n", length(timepoints)))
    cat(sprintf("  True parents: %d\n", length(parent_barcodes)))
    cat("================================================================\n\n")
  }
  
  all_clusters <- list()
  all_validations <- list()
  summary_df <- data.frame()
  
  for (tp in timepoints) {
    
    barcode_data <- sim_result$timepoint_data[[tp]]
    
    cluster_input <- barcode_data %>%
      dplyr::select(barcode, counts)
    
    cluster_result <- super_cluster2(
      cluster_input,
      distance = distance,
      method = method,
      verbose = FALSE
    )
    
    validation <- evaluate_single_clustering(
      cluster_result,
      barcode_data,
      parent_barcodes
    )
    
    all_clusters[[tp]] <- cluster_result
    all_validations[[tp]] <- validation
    
    tp_num <- barcode_data$timepoint_num[1]
    
    summary_row <- data.frame(
      timepoint = tp,
      timepoint_num = tp_num,
      n_barcodes = nrow(barcode_data),
      n_clusters = nrow(cluster_result),
      n_parents = length(parent_barcodes),
      purity = validation$purity,
      completeness = validation$completeness,
      f1_score = validation$f1_score,
      parent_recovery = validation$parent_recovery_rate
    )
    
    summary_df <- dplyr::bind_rows(summary_df, summary_row)
    
    if (verbose) {
      cat(sprintf("%s: %3d barcodes -> %2d clusters | Purity: %.2f | Complete: %.2f\n",
                  tp, nrow(barcode_data), nrow(cluster_result),
                  validation$purity, validation$completeness))
    }
  }
  
  if (verbose) {
    cat("\n")
    cat("----------------------------------------------------------------\n")
    cat(sprintf("OVERALL AVERAGE:  Purity: %.3f | Completeness: %.3f | F1: %.3f\n",
                mean(summary_df$purity),
                mean(summary_df$completeness),
                mean(summary_df$f1_score)))
    cat("----------------------------------------------------------------\n")
  }
  
  return(list(
    clusters = all_clusters,
    validations = all_validations,
    summary = summary_df,
    distance = distance,
    method = method,
    sim_result = sim_result
  ))
}

# -----------------------------------------------------------------------------
# Compare distance thresholds
# -----------------------------------------------------------------------------

#' @export
compare_distance_thresholds <- function(sim_result,
                                        distances = c(1, 2, 3, 4, 5),
                                        verbose = TRUE) {
  
  if (verbose) {
    cat("\n")
    cat("================================================================\n")
    cat("     COMPARING DISTANCE THRESHOLDS                             \n")
    cat("================================================================\n\n")
  }
  
  all_results <- list()
  comparison_df <- data.frame()
  
  for (d in distances) {
    if (verbose) cat(sprintf("Testing distance = %d...\n", d))
    
    result <- cluster_and_validate_timeseries(sim_result, distance = d, verbose = FALSE)
    
    all_results[[as.character(d)]] <- result
    
    summary_row <- data.frame(
      distance = d,
      mean_purity = mean(result$summary$purity),
      mean_completeness = mean(result$summary$completeness),
      mean_f1 = mean(result$summary$f1_score),
      mean_parent_recovery = mean(result$summary$parent_recovery),
      mean_n_clusters = mean(result$summary$n_clusters)
    )
    
    comparison_df <- dplyr::bind_rows(comparison_df, summary_row)
  }
  
  if (verbose) {
    cat("\n")
    cat("Distance | Purity | Complete | F1     | Recovery | Clusters\n")
    cat("---------|--------|----------|--------|----------|----------\n")
    for (i in 1:nrow(comparison_df)) {
      cat(sprintf("   %d     |  %.3f |   %.3f  |  %.3f |   %.3f  |   %.1f\n",
                  comparison_df$distance[i],
                  comparison_df$mean_purity[i],
                  comparison_df$mean_completeness[i],
                  comparison_df$mean_f1[i],
                  comparison_df$mean_parent_recovery[i],
                  comparison_df$mean_n_clusters[i]))
    }
    
    best_f1_idx <- which.max(comparison_df$mean_f1)
    cat(sprintf("\n-> Best F1 score at distance = %d\n", comparison_df$distance[best_f1_idx]))
  }
  
  return(list(results = all_results, comparison = comparison_df))
}

# -----------------------------------------------------------------------------
# Plotting functions
# -----------------------------------------------------------------------------

#' @export
plot_validation_metrics <- function(validation_result) {
  
  summary_df <- validation_result$summary
  
  metrics_long <- summary_df %>%
    dplyr::select(timepoint_num, purity, completeness, f1_score, parent_recovery) %>%
    tidyr::pivot_longer(
      cols = c(purity, completeness, f1_score, parent_recovery),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(metric = factor(metric, 
                                  levels = c("purity", "completeness", "f1_score", "parent_recovery"),
                                  labels = c("Purity", "Completeness", "F1 Score", "Parent Recovery")))
  
  ggplot2::ggplot(metrics_long, ggplot2::aes(x = timepoint_num, y = value, color = metric)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      title = sprintf("Clustering Validation (distance = %d)", validation_result$distance),
      x = "Time Point", y = "Score", color = "Metric"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' @export
plot_distance_comparison <- function(comparison_result) {
  
  df <- comparison_result$comparison %>%
    tidyr::pivot_longer(
      cols = c(mean_purity, mean_completeness, mean_f1, mean_parent_recovery),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(metric = factor(metric,
                                  levels = c("mean_purity", "mean_completeness", "mean_f1", "mean_parent_recovery"),
                                  labels = c("Purity", "Completeness", "F1 Score", "Parent Recovery")))
  
  ggplot2::ggplot(df, ggplot2::aes(x = distance, y = value, color = metric)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::scale_x_continuous(breaks = unique(df$distance)) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      title = "Clustering Performance vs Distance Threshold",
      x = "Distance Threshold", y = "Score", color = "Metric"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' @export
plot_true_vs_clustered <- function(validation_result) {
  
  sim_result <- validation_result$sim_result
  
  # True abundances
  true_abundance <- as.data.frame(sim_result$abundance_proportions)
  colnames(true_abundance) <- sim_result$parent_barcodes
  true_abundance$timepoint_num <- 0:(nrow(true_abundance) - 1)
  
  true_long <- true_abundance %>%
    tidyr::pivot_longer(cols = -timepoint_num, names_to = "parent", values_to = "proportion") %>%
    dplyr::mutate(source = "True Abundances")
  
  parent_order <- true_long %>%
    dplyr::group_by(parent) %>%
    dplyr::summarise(total = sum(proportion), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total)) %>%
    dplyr::pull(parent)
  
  true_long <- true_long %>%
    dplyr::mutate(parent = factor(parent, levels = rev(parent_order)))
  
  # Clustered abundances
  clusters_list <- validation_result$clusters
  validations_list <- validation_result$validations
  
  clustered_data <- data.frame()
  for (tp in names(clusters_list)) {
    cluster_result <- clusters_list[[tp]]
    validation <- validations_list[[tp]]
    tp_num <- as.numeric(gsub("T", "", tp))
    
    for (i in 1:nrow(cluster_result)) {
      clustered_data <- dplyr::bind_rows(clustered_data, data.frame(
        timepoint_num = tp_num,
        parent = validation$cluster_dominant_parent[i],
        counts = cluster_result$sum_counts[i]
      ))
    }
  }
  
  clustered_long <- clustered_data %>%
    dplyr::group_by(timepoint_num, parent) %>%
    dplyr::summarise(counts = sum(counts), .groups = "drop") %>%
    dplyr::group_by(timepoint_num) %>%
    dplyr::mutate(proportion = counts / sum(counts)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(source = "Clustered") %>%
    dplyr::select(-counts) %>%
    dplyr::mutate(parent = factor(parent, levels = levels(true_long$parent)))
  
  # Combine
  combined <- dplyr::bind_rows(true_long, clustered_long) %>%
    dplyr::mutate(source = factor(source, levels = c("True Abundances", "Clustered")))
  
  n_parents <- length(parent_order)
  colors <- RColorBrewer::brewer.pal(max(3, min(n_parents, 12)), "Set3")
  if (n_parents > 12) colors <- colorRampPalette(colors)(n_parents)
  
  ggplot2::ggplot(combined, ggplot2::aes(x = timepoint_num, y = proportion, fill = parent)) +
    ggplot2::geom_area(alpha = 0.8, color = "white", linewidth = 0.2) +
    ggplot2::facet_wrap(~source, ncol = 1) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = "True vs Clustered Abundances",
      x = "Time Point", y = "Relative Abundance", fill = "Parent"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 6, family = "mono")
    )
}