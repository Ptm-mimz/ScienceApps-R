```{r}
library(cellpaintr)
library(ggplot2)
library(scater)
library(scran)
library(bluster)
library(tidyr)
library(dplyr)
```

# path to files

```{r}

file_plate <- paste0(
  "/scratch/plapha/250320 stress panel_cytokine/output/",
  "MyExpt_Cells.csv"
)

sce <- loadData(file_plate)

```

## Data Transformation
```{r}
plotCellsPerImage(sce)
sce <- removeOutliers(sce, min = 10, max = 300)

#sce <- removeMissingValues(sce)
# remove missing values in cells
mat <- assay(sce, name = "features")
miss_perc_cells <- apply(mat, 2, function(x) mean(is.na(x)))
cell_ids <- which(miss_perc_cells > 0)
sce <- sce[, -cell_ids]

sce <- removeLowVariance(sce)
sce <- transformLogScale(sce)

plotCellsPerImage(sce)
sce <- removeOutliers(sce, min = 10, max = 300)
 
 
sce <- removeLowVariance(sce)
sce <- transformLogScale(sce)
```

```{r}
# sce1 <- sce
sce <- sce1

#if I want to delete feature that has 0 value
mat <- assay(sce, name = "features")
miss_perc_features <- apply(mat, 1, function(x) mean(is.na(x)))
feature_ids <- which(miss_perc_features > 0)
sce <- sce[-feature_ids,]

###############################################################################


#if I want to delete cell that has 0 value, but keep feature data
mat <- assay(sce, name = "features")
miss_perc_cells <- apply(mat, 2, function(x) mean(is.na(x)))
cell_ids <- which(miss_perc_cells > 0)
sce <- sce[, -cell_ids]

###############################################################################

sce <- removeLowVariance(sce)
sce <- transformLogScale(sce)

```


filter only full stained

```{r}

# Define the single-stained wells
single_stained_wells <- c("B05", "C05", "D05", "E05", "F05")

# Extract the well information from colData
well_info <- colData(sce)$Well

# Filter out single-stained wells
full_stained_indices <- !(well_info %in% single_stained_wells)

# Subset the sce object to retain only full-stained data
sce_full_stained <- sce[, full_stained_indices]

````
filter only Intensity data

````{r}
feature_names <- rownames(sce)

#see all intensity feature name
intensity_features <- grep("^Intensity", feature_names, value = TRUE)
print(intensity_features)


# Extract intensity features from the assay slot
intensity_data <- as.data.frame(t(assay(sce_full_stained, "features")))  # change back to "tfmfeatures" if work with logtransform in above step

# Filter for columns that start with "Intensity"
intensity_features <- grep("^Intensity", colnames(intensity_data), value = TRUE)
#intensity_data <- intensity_data[, intensity_features, drop = FALSE] #if want all Intensity data

intensity_data <- intensity_data[, c("Intensity_MeanIntensity_CorrLyso", 
                                    "Intensity_MeanIntensity_CorrTMRM", 
                                    "Intensity_MeanIntensity_CorrCytoID",
                                    "Intensity_IntegratedIntensity_CorrLyso",
                                    "Intensity_IntegratedIntensity_CorrTMRM",
                                    "Intensity_IntegratedIntensity_CorrCytoID" )]


# Extract metadata from colData
metadata <- as.data.frame(colData(sce_full_stained)[, c("Patient", "Treatment", "Disease")])  # Replace with correct column names

# Combine metadata with intensity data
df <- cbind(metadata, intensity_data)

df_long <- pivot_longer(
  df,
  cols = starts_with("Intensity"),
  names_to = "Feature",
  values_to = "Intensity"
)

# Aggregate intensity values by Patient and Treatment
df_aggregated <- df_long %>%
  group_by(Patient, Treatment, Feature, Disease) %>%
  summarize(Intensity = mean(Intensity, na.rm = TRUE)) %>%  # Use mean or median
  ungroup()

````

# Create the box plot

````{r}
# Create the plot with Healthy patients in blue
ggplot(df_aggregated, aes(x = Treatment, y = Intensity, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +  # Filled boxes with transparency
  geom_jitter(aes(color = Disease), width = 0.2, alpha = 0.5, size = 2) +  # Map 'Disease' to color
  scale_color_manual(values = c("Healthy" = "blue", "RA" = "red")) +  # Set colors for 'Disease'
  facet_wrap(~Feature, scales = "free_y") +  # Facet by feature
  labs(
    title = "Comparison of Intensity Features Across Treatments",
    subtitle = "Each dot represents one patient (Healthy in blue)",
    x = "Treatment",
    y = "Intensity",
    fill = "Treatment",  # Legend for fill (boxplot)
    color = "Disease"    # Legend for color (dots)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 14),  # Increase general text size
    axis.text = element_text(size = 12),  # Increase axis text size
    axis.title = element_text(size = 14),  # Increase axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Increase plot title size and make it bold
    plot.subtitle = element_text(size = 14),  # Increase plot subtitle size
    strip.text = element_text(size = 14)  # Increase facet label size
  )

````
Add stat comparison

````{r}
library(ggplot2)
library(ggpubr)
library(multcompView)

# Run ANOVA
anova_result <- aov(Intensity ~ Treatment, data = df_aggregated)

# Run Tukey's HSD test
tukey_result <- TukeyHSD(anova_result)

# Extract significant comparisons from Tukey's HSD
tukey_df <- as.data.frame(tukey_result$Treatment)
tukey_df$Comparison <- rownames(tukey_df)
tukey_df$Significance <- ifelse(tukey_df$`p adj` < 0.05, "*", "ns")

# Manually extract Treatment groups
tukey_df$Group1 <- sapply(strsplit(tukey_df$Comparison, "-"), `[`, 1)
tukey_df$Group2 <- sapply(strsplit(tukey_df$Comparison, "-"), `[`, 2)

# Create the plot with ANOVA and Tukey's HSD
ggplot(df_aggregated, aes(x = Treatment, y = Intensity, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +  # Filled boxes with transparency
  geom_jitter(width = 0.2, alpha = 0.5, size = 2, color = "black") +  # Add dots for each patient
  facet_wrap(~Feature, scales = "free_y") +  # Facet by feature
  stat_compare_means(method = "anova", label.y = max(df_aggregated$Intensity, na.rm = TRUE) * 1.1) +  # Add ANOVA p-value
  labs(
    title = "Comparison of Intensity Features Across Treatments",
    subtitle = "Each dot represents one patient",
    x = "Treatment",
    y = "Intensity",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text = element_text(size = 14)
  ) +
  geom_text(data = tukey_df, 
            aes(x = (as.numeric(factor(Group1)) + as.numeric(factor(Group2))) / 2, 
                y = max(df_aggregated$Intensity, na.rm = TRUE) * 1.15, 
                label = Significance), 
            inherit.aes = FALSE, 
            size = 5)

````
## Unsupervised Analysis

Use `scater` for exploratory data analysis.

PCA:
````{r}
sce <- runPCA(sce, exprs_values = "tfmfeatures")
plotReducedDim(sce, dimred = "PCA", colour_by = "Disease")
plotReducedDim(sce, dimred = "PCA", colour_by = "Treatment")

```
Custom plots using `ggcells` function.

```{r}
plotReducedDim(sce, dimred = "PCA", colour_by = "Treatment")
    ggcells(sce, aes(x = PCA.2, y = PCA.3, color = Treatment)) +
      geom_point() +
      theme_minimal() +
      labs(title = "PCA: PC2 vs PC3", x = "PC2", y = "PC3")
```

UMAP:

```{r}
sce <- runUMAP(sce, exprs_values = "tfmfeatures")
plotReducedDim(sce, dimred = "UMAP", colour_by = "Disease")
plotReducedDim(sce, dimred = "UMAP", colour_by = "Treatment")
```

Faster clustering:

```{r}
set.seed(23)
kgraph.clusters <- clusterCells(sce_full_stained, use.dimred = "PCA",
    BLUSPARAM = TwoStepParam(
        first = KmeansParam(centers = 2),
        second = SNNGraphParam(k = 2, type = "rank", cluster.fun = "walktrap")
    )
)
table(kgraph.clusters)
colLabels(sce_full_stained) <- kgraph.clusters
plotReducedDim(sce_full_stained, "UMAP", colour_by = "label")
```


## Fit Multiple Models

Train machine learning models on multiple perturbations.

```{r}
interest <- setdiff(unique(sce$Treatment), c("cont"))
result_list <- lapply(interest, function(level) {
  fitModel(sce, 
           target = "Treatment", 
           interest_level = level, reference_level = "cont", 
           group = "Patient", strata = c("Disease"), 
           n_folds = 20, n_threads = 8)
})
```

Plot AUC comparison between perturbations.

```{r}
plotAUC(result_list)
```

