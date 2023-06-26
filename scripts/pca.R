library(rsvd) # install.packages("rsvd")
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")

# Replace the file paths and delimiters as needed
file_path <- "/lizardfs/guarracino/keane_mouse_pangenome/short_reads/alignments/matrix.50M.tsv.gz"
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

# Write the data frame to a CSV file
write.csv(principal_components_df, paste0(file_path, ".pc_1to2.csv"), row.names = FALSE)
# When you want to read the data back into R:
#principal_components_df <- read.csv(paste0(file_path, ".pc_1to2.csv"))

# Plot the first two principal components using ggplot2
p <- ggplot(principal_components_df, aes(x = X1, y = X2, label = Sample)) +
  geom_point() +
  geom_text_repel(size = 3, max.overlaps = 1000) + # Use geom_text_repel instead of geom_text, and adjust the size if needed
  labs(x = paste0("PC1: ", round(explained_variance_ratio[1] * 100, 2), "% variance"), 
       y = paste0("PC2: ", round(explained_variance_ratio[2] * 100, 2), "% variance")) +
  theme_bw() +
  expand_limits(x = c(min(principal_components_df$X1), max(principal_components_df$X1)), 
                y = c(min(principal_components_df$X2), max(principal_components_df$X2))) # Use expand_limits to ensure all labels fit in the plot

ggsave(paste0(file_path, ".pc_1to2.pdf"), width = 20, height = 12, plot = p)
