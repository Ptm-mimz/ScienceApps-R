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
# 
# file_plate <- paste0(
#   "/data/plapha/Cellprofiler analysis result files/250228 lipid panel_hand/",
#   "MyExpt_Cells.csv"
# )

 file_plate1 <- paste0(
   "/data/plapha/Cellprofiler analysis result files/250709 cytokines_metabolism/",
   "MyExpt_Cells-plate1.csv"
 )

 file_plate2 <- paste0(
   "/data/plapha/Cellprofiler analysis result files/250709 cytokines_metabolism/",
   "MyExpt_Cells-plate2.csv"
 )

  file_plate3 <- paste0(
   "/data/plapha/Cellprofiler analysis result files/250717 cytokines_metabolism lot2/",
   "MyExpt_Cells_plate1.csv"
 )
  
 
  file_plate4 <- paste0(
   "/data/plapha/Cellprofiler analysis result files/250717 cytokines_metabolism lot2/",
   "MyExpt_Cells_plate2.csv"
 )
  
  sce1 <- loadData(file_plate1)
  sce2 <- loadData(file_plate2)
  sce3 <- loadData(file_plate3)
  sce4 <- loadData(file_plate4)
  
  sce <- cbind(sce1, sce2, sce3, sce4)
  
```
Prepare the data for machine learning.

```{r}
# fix plate id
#sce$Plate <- unique(sce$Patient)[1]
 
# remove missing values in cells
mat <- assay(sce, name = "features")
miss_cells <- apply(mat, 2, function(x) sum(is.na(x)))
cell_ids <- which(miss_cells > 0)
#sce <- sce[, -cell_ids] #skip becases there was no missing value

 
# define reference level
sce$Treatment <- as.factor(sce$Treatment)
sce$Treatment <- relevel(sce$Treatment, ref = "control")
 
# batch correct and normalize features
sce <- normalizeExclude(sce)
 
# remove outliers
plotCellsPerImage(sce)
sce <- removeOutliers(sce, min = 50, max = 300)

```

## Fit Single Model
```{r}
result <- fitModel(sce,
                   target = "Disease", 
                   interest_level = "RA", reference_level = "healthy", 
                   group = "Patient", strata =  "Disease", 
                   n_folds = 20, n_threads = 8)

plotROC(result)


```

For Multiple model

```{r}

interest <- setdiff(unique(sce$Treatment), c("control"))
result_list <- lapply(interest, function(level) {
  fitModel(sce, 
           target = "Treatment", 
           interest_level = level, reference_level = "control", 
           group = "Patient", strata = c("Disease"), 
           n_folds = 10, n_threads = 8)
})

plotAUC(result_list)


```
Only selected features

````{r}

# overview of CellProfiler feature categories
#   include: Texture, Intensity, AreaShape, 
#   maybe include: RadialDistribution,  Neighbors, Location
#   exclude because many corner cases: Granularity, Correlation

sce_select <- sce[str_detect(names(sce), "Texture|Intensity|AreaShape"), ]

#sce_select <- sce_select[str_detect(names(sce_select), "CorrCyto"), ]

result <- fitModel(sce_select,
                   target = "Disease", 
                   interest_level = "RA", reference_level = "healthy", 
                   group = "Patient", strata =  "Disease", 
                   n_folds = 20, n_threads = 8)

plotROC(result)


```

filter only Intensity data

````{r}
feature_names <- rownames(sce)

#see all intensity feature name
intensity_features <- grep("Intensity", feature_names, value = TRUE)
print(intensity_features)


# Extract intensity features from the assay slot
intensity_data <- as.data.frame(t(assay(sce, "features")))  # change back to "tfmfeatures" if work with logtransform in above step

# Filter for columns that start with "Intensity"
intensity_features <- grep("^Intensity", colnames(intensity_data), value = TRUE)
#intensity_data <- intensity_data[, intensity_features, drop = FALSE] #if want all Intensity data

intensity_data <- intensity_data[, c("Intensity_MeanIntensity_CorrMitosox", 
                                    "Intensity_MeanIntensity_CorrLipid", 
                                    "Intensity_MeanIntensity_CorrEdU", 
                                    "Intensity_IntegratedIntensity_CorrMitosox",
                                    "Intensity_IntegratedIntensity_CorrLipid",
                                    "Intensity_IntegratedIntensity_CorrEdU" )]


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

Handle with EdU plot
````{r}
# Plot histogram for Intensity_IntegratedIntensity_CorrEdU

# Define the cutoff value
cutoff <- 0.001

# Create the histogram with adjusted axes
ggplot(df, aes(x = Intensity_MeanIntensity_CorrEdU)) +
  geom_histogram(binwidth = 0.0001, fill = "blue", alpha = 0.7) +  # Adjust binwidth as needed
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1) +  # Add cutoff line
  labs(
    title = "Histogram of Intensity_MeanIntensity_CorrEdU",
    x = "Intensity_MeanIntensity_CorrEdU",
    y = "Count"
  ) +
  facet_wrap(Disease ~ Treatment) +  # Facet by Disease and Treatment
  scale_x_continuous(
    limits = c(0, 0.005),  # Set x-axis limits from 0 to 0.005
    breaks = seq(0, 0.005, by = 0.001)  # Set x-axis breaks every 0.001
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),  # Increase general text size
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Increase x-axis label size and rotate
    axis.text.y = element_text(size = 12),  # Increase y-axis label size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Increase plot title size and make it bold
    strip.text = element_text(size = 14)  # Increase facet label size
  )

result <- df %>%
  group_by(Treatment, Disease) %>%
  summarize(
    above_cutoff_mean = sum(Intensity_MeanIntensity_CorrEdU > cutoff),  # Number of cells above cutoff
    percentage_above_cutoff_mean = (above_cutoff_mean / n()) * 100  # Percentage of cells above cutoff
  ) %>%
  ungroup()

# Print the result
print(result)

#plot percentage on histogram
ggplot(df, aes(x = Intensity_MeanIntensity_CorrEdU)) +
  geom_histogram(binwidth = 0.0001, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Histogram of Intensity_MeanIntensity_CorrEdU",
    x = "Intensity_MeanIntensity_CorrEdU",
    y = "Count"
  ) +
  facet_wrap(Disease ~ Treatment) +
  scale_x_continuous(
    limits = c(0, 0.005),  # Set x-axis limits from 0 to 0.005
    breaks = seq(0, 0.005, by = 0.001)  # Set x-axis breaks every 0.001
  ) +
  geom_text(
    data = result,  # Use the result data frame for annotations
    aes(x = 0.004, y = Inf, label = paste("Above cutoff:", above_cutoff_mean, "\n(", round(percentage_above_cutoff_mean, 2), "%)")),
    vjust = 2, hjust = 1, color = "black", size = 4
  ) +
  theme_minimal()
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
    subtitle = "Each dot represents one patient (RA is in red)",
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
Plot patient box side-by-side

````{r}
# Create the plot with Healthy and RA patients side-by-side
ggplot(df_long, aes(x = Treatment, y = Intensity, fill = Treatment)) +
  geom_boxplot(aes(color = Disease), outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.8)) +  # Boxplots with side-by-side grouping
  geom_jitter(aes(color = Disease), alpha = 0.5, size = 2, position = position_dodge(width = 0.8)) +  # Jittered points with side-by-side grouping (no width argument)
  scale_color_manual(values = c("Healthy" = "blue", "RA" = "red")) +  # Set colors for 'Disease'
  scale_fill_manual(values = c("Treatment1" = "lightgray", "Treatment2" = "darkgray")) +  # Set fill colors for treatments
  facet_wrap(~Feature, scales = "free_y") +  # Facet by feature
  labs(
    title = "Comparison of Intensity Features Across Treatments",
    subtitle = "Each dot represents one patient (Healthy in blue, RA in red)",
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
Plot hist of Lipidtox to set threshold

````{r}
#### have to go back to load the raw sheet again(before log transform)
# Define the cutoff value

cutoff <- 1000

# Create the histogram with adjusted axes
ggplot(df, aes(x = Intensity_IntegratedIntensity_CorrLipid)) +
  geom_histogram(binwidth = 5, fill = "blue", alpha = 0.7) +  # Adjust binwidth as needed
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", linewidth = 1) +  # Add cutoff line
  labs(
    title = "Histogram of Intensity_IntegratedIntensity_CorrLipid",
    x = "Intensity_IntegratedIntensity_CorrLipid",
    y = "Count"
  ) +
  facet_wrap(Disease ~ Treatment) +  # Facet by Disease and Treatment
  scale_x_continuous(
    limits = c(0, 10000),  # Set x-axis limits from 0 to 0.005
    breaks = seq(0, 10000, by = 1000)  # Set x-axis breaks every 0.001 if data are smaller number
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),  # Increase general text size
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Increase x-axis label size and rotate
    axis.text.y = element_text(size = 12),  # Increase y-axis label size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Increase plot title size and make it bold
    strip.text = element_text(size = 14)  # Increase facet label size
  )

#calculate percentage
result <- df %>%
  group_by(Treatment, Disease) %>%
  summarize(
    above_cutoff_mean = sum(Intensity_IntegratedIntensity_CorrLipid > cutoff),  # Number of cells above cutoff
    percentage_above_cutoff_mean = (above_cutoff_mean / n()) * 100  # Percentage of cells above cutoff
  ) %>%
  ungroup()

# Print the result
print(result)

#Show percentage on histogram
ggplot(df, aes(x = Intensity_IntegratedIntensity_CorrLipid)) +
  geom_histogram(binwidth = 5, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Histogram of Intensity_IntegratedIntensity_CorrLipid",
    x = "Intensity_IntegratedIntensity_CorrLipid",
    y = "Count"
  ) +
  facet_wrap(Disease ~ Treatment) +
  scale_x_continuous(
    limits = c(0, 10000),  # Set x-axis limits from 0 to 0.005
    breaks = seq(0, 10000, by = 1000)  
  ) +
  geom_text(
    data = result,  # Use the result data frame for annotations
    aes(x = 8000, y = Inf, label = paste("Above cutoff:", above_cutoff_mean, "\n(", round(percentage_above_cutoff_mean, 2), "%)")),
    vjust = 2, hjust = 1, color = "black", size = 4
  ) +
  theme_minimal()



````
## Fit Multiple Models

Train machine learning models on multiple perturbations.

```{r}
interest <- setdiff(unique(sce_full_stained$Treatment), c("cont"))
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

````{r}



````
