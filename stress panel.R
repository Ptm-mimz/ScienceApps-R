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

Prepare the data for machine learning.
```{r}
# fix plate id
sce$Plate <- unique(sce$Patient)[1]
 
# remove missing values in cells
mat <- assay(sce, name = "features")
miss_cells <- apply(mat, 2, function(x) sum(is.na(x)))
cell_ids <- which(miss_cells > 0)
sce <- sce[, -cell_ids]
 
# define reference level
sce$Treatment <- as.factor(sce$Treatment)
sce$Treatment <- relevel(sce$Treatment, ref = "control")
 
# batch correct and normalize features
sce <- normalizeExclude(sce)
 
# remove outliers
plotCellsPerImage(sce)
sce <- removeOutliers(sce, min = 20, max = 500)

```
flter Intensity data

```{r}
#EXTRACT NORMALISED FEATURES
mat <- assay(sce, "features")

#Get only 'Intensity' features
intensity_features <- grep("^Intensity_", rownames(mat), value = TRUE)
mat_intensity <- mat[intensity_features, ]

#transpose for PCA
mat_pca <- t(mat_intensity)

#run PCA
pca <- prcomp(mat_pca, scale. = TRUE)

# Calculate % variance for each PC
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_label <- paste0("PC1 (", round(percentVar[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(percentVar[2], 1), "%)")

#plot PCA
pca_df <- as.data.frame(pca$x)
pca_df$Treatment <- sce$Treatment  # add metadata if needed

ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(alpha = 0.3, size = 1.5) +
  theme_minimal() +
  labs(
    title = "PCA of Intensity Features",
    x = pc1_label,
    y = pc2_label
  ) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

```
UMAp from Intensity data

```{r}

library(uwot)

mat_umap <- t(mat_intensity)

set.seed(42)  # for reproducibility
umap_result <- umap(mat_umap, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")

#plot UMAP
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Treatment <- sce$Treatment


ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Treatment)) +
  geom_point(alpha = 0.3, size = 1.5) +
  facet_wrap(~ Treatment) +
  theme_minimal() +
  labs(
    title = "UMAP of Intensity Features",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

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
intensity_data <- as.data.frame(t(assay(sce, "tfmfeatures")))  # change back to "tfmfeatures" if work with logtransform in above step

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
metadata <- as.data.frame(colData(sce)[, c("Patient", "Treatment", "Disease")])  # Replace with correct column names

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

```{r}

sce <- runPCA(sce, exprs_values = "tfmfeatures")
plotReducedDim(sce, dimred = "PCA", colour_by = "Patient")
plotReducedDim(sce, dimred = "PCA", colour_by = "Treatment")
plotPCACor(sce, filter_by = 1, top = 20)
plotPCACor(sce, filter_by = 2, top = 20)
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
plotReducedDim(sce, dimred = "UMAP", colour_by = "Patient")
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
interest <- setdiff(unique(sce$Treatment), c("control"))
sce$Treatment <- as.factor(sce$Treatment)
sce$Treatment <- droplevels(sce$Treatment)

result_list <- lapply(interest, function(level) {
  fitModel(sce,
           target = "Treatment",
           interest_level = level, reference_level = "control",
           group = "Patient", strata = c("Disease"),
           n_folds = 20, n_threads = 8)
})

plotAUC(result_list)


##test diff if panel can distinguish between disease

sce_cont <- sce[, sce$Treatment == "control"]
 
result <- fitModel(sce_cont,
                   target = "Disease", 
                   interest_level = "RA", reference_level = "Healthy", 
                   group = "Patient", strata = "Disease",
                   n_folds = 20, n_threads = 8)



#in healthy with treatment _ how does it look?
sce_healthy <- sce[, sce$Disease == "healthy"]
 
result <- fitModel(sce_healthy,
                   target = "Treatment", 
                   interest_level = "PIC", reference_level = "control", 
                   group = "Patient", strata = "Disease",
                   n_folds = 20, n_threads = 8)
```

Plot AUC comparison between perturbations.

```{r}
plotROC(result)

```
Check sample info of individual curves.

```{r}
# identify fold id
plotROC(result) + 
  facet_wrap(~fold) + # plot each curves on separate facet
  geom_path(alpha = 1.0) # for better visibility

# print sample info for one fold
printROC(result, 11)

```
