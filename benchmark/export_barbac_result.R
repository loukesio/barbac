library(barbac)
library(dplyr)
library(tidyr)
library(readr)

work_dir  <- "~/Documents/Projects/Barcodes/barbac-benchmark"
input_csv <- file.path(work_dir, "barbac_benchmark_input.csv")

result <- super_cluster2(
  input_csv,
  distance    = 3,
  merge_ratio = 20
)

centroids <- data.frame(
  central_barcode = result$central_barcode,
  sum_counts      = result$sum_counts
)
write_csv(centroids, file.path(work_dir, "barbac_result.csv"))
cat("Wrote", nrow(centroids), "centroids to barbac_result.csv\n")

members <- result %>%
  select(central_barcode, all_barcodes, all_counts) %>%
  mutate(row = row_number()) %>%
  unnest(c(all_barcodes, all_counts)) %>%
  transmute(
    central_barcode,
    member       = all_barcodes,
    member_count = all_counts
  )
write_csv(members, file.path(work_dir, "barbac_clusters.csv"))
cat("Wrote", nrow(members), "members to barbac_clusters.csv\n")
