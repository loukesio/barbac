library(dplyr)

sample_table <- tibble::tibble(
  sample = c("sample1", "sample2"),
  R1 = c("path/to/sample1_R1.fastq", "path/to/sample2_R1.fastq"),
  R2 = c("path/to/sample1_R2.fastq", "path/to/sample2_R2.fastq")
)


sample_table %>% 
  write.csv(., "example_sample_table.csv", row.names = FALSE)
