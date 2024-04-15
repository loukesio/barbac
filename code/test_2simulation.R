data <- data.frame(
  error_rate = c(0, 0.1, 0.2, 0.5),
  average_clustering_success_percentage = c(100, 99.8, 99.6, 79.3),
  standard_deviation = c(0, 0.134, 0.1589, 0.17853)
)

library(lukesky)

data %>% 
  ggplot() +
  geom_col(aes(x=error_rate, average_clustering_success_percentage), fill="#F6AB5A") +
   theme_bw() +
  theme_volcano() +
  labs(x="Error rate", y="Average clustering success",
        title="Average clustering success after testing\nbarbac_cluster on 100 simualted datasets")
