# Load required packages
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(FactoMineR)
library(factoextra)

# Read the data (assumes CSV, adjust if needed)
df <- read.csv("your_data.csv", stringsAsFactors = FALSE)

# 1. Histogram of Best Localization Probability
ggplot(df, aes(x = Best.Localization.Probability)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  theme_minimal() +
  xlab("Best Localization Probability") +
  ylab("Number of Peptides") +
  ggtitle("Distribution of Localization Probability")

# 2. Normalize intensity values (log2 transform + pseudocount)
intensity_cols <- grep("Intensity", names(df), value = TRUE)
df[intensity_cols] <- log2(df[intensity_cols] + 1)

# 3. Heatmap of intensities
# Subset matrix for heatmap
heatmap_mat <- df %>%
  select(all_of(intensity_cols)) %>%
  as.matrix()
rownames(heatmap_mat) <- df$Index

# Optional: remove rows with all 0s or low variance
heatmap_mat <- heatmap_mat[rowSums(heatmap_mat) > 0, ]

# Plot heatmap
pheatmap(heatmap_mat,
         scale = "row",
         show_rownames = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Peptide Expression Heatmap")

# 4. PCA
pca_data <- t(heatmap_mat)
pca_result <- prcomp(pca_data, scale. = TRUE)
fviz_pca_ind(pca_result,
             repel = TRUE,
             title = "PCA of Samples Based on Peptide Intensities")

# 5. Bar plot of most frequent localizations (if categorical labels exist)
# If you have categorical localization assignments, plot them
# Example assumes "Best Localization Probability" is already grouped (e.g., "Nucleus", "Cytoplasm")
# If not, bin it manually or add category column

# 6. Volcano plot (only if you have p-values and fold changes)
# Not available in your table directly; would need statistical testing first

# Save plots if needed
ggsave("localization_histogram.png")