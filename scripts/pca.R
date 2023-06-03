library(rsvd) # install.packages("rsvd")
library(tidyverse) # install.packages("tidyverse")

# Replace the file paths and delimiters as needed
file_path <- "/lizardfs/guarracino/keane_mouse_pangenome/short_reads/alignments/matrix.tsv.gz"
sample_names_file <- "/lizardfs/guarracino/keane_mouse_pangenome/short_reads/alignments/samples.txt"

data_matrix <- read.csv(file_path, sep = "\t", header = FALSE)
sample_names <- read.csv(sample_names_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

## Remove rows with all zeros
#data_matrix <- data_matrix[rowSums(data_matrix != 0, na.rm = TRUE) > 0, ]

data_matrix <- scale(t(data_matrix))

randomized_svd_result <- rsvd(data_matrix, k = 2)

# Compute principal components
principal_components <- data_matrix %*% randomized_svd_result$v

# Compute variance explained by each principal component
explained_variance_ratio <- randomized_svd_result$d^2 / sum(randomized_svd_result$d^2)

# Convert principal components to a data frame and add sample names
principal_components_df <- data.frame(principal_components[, 1:2]) %>%
  mutate(Sample = sample_names$V1)

# Plot the first two principal components using ggplot2
p <- ggplot(principal_components_df, aes(x = X1, y = X2, label = Sample)) +
  geom_point() +
  geom_text(vjust = 1, hjust = 1, size = 3) +
  labs(x = paste0("PC1: ", round(explained_variance_ratio[1] * 100, 2), "% variance"), 
       y = paste0("PC2: ", round(explained_variance_ratio[2] * 100, 2), "% variance")) +
  theme_bw()

ggsave("pca.pdf", plot = p)
