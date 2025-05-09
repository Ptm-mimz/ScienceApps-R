```{r}
library(cellpaintr)
library(ggplot2)
library(scater)
library(scran)
library(bluster)
library(tidyr)
library(dplyr)
library(stringr)
```

# path to files

```{r}

file_plate1 <- paste0(
  "/scratch/plapha/CP_241211 cellpainting shoulder/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate2 <- paste0(
  "/scratch/plapha/CP_241219 cellpainting hands/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate3 <- paste0(
  "/scratch/plapha/CP_240828 shoulder all drugtx/plate1/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate4 <- paste0(
  "/scratch/plapha/CP_240828 shoulder all drugtx/plate2/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate5 <- paste0(
  "/scratch/plapha/CP_240828 shoulder all drugtx/plate3/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate6 <- paste0(
  "/scratch/plapha/CP_240828 shoulder all drugtx/plate4/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate7 <- paste0(
  "/scratch/plapha/CP_240925 shoulder lot2/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate8 <- paste0(
  "/scratch/plapha/CP_240820 OA knee/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate9 <- paste0(
  "/scratch/plapha/CP_240805 knee plate1 reimage/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate10 <- paste0(
  "/scratch/plapha/CP_240805 knee plate2/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate11 <- paste0(
  "/scratch/plapha/CP_240805 shoulder batch2/output-newplatename/",
  "MyExpt_Cells_2.csv"
)

file_plate12 <- paste0(
  "/scratch/plapha/CP_240802 shoulder batch1/plate1/output-newplatename2/",
  "MyExpt_Cells.csv"
)

file_plate13 <- paste0(
  "/scratch/plapha/CP_240802 shoulder batch1/plate2/output-newplatename/",
  "MyExpt_Cells.csv"
)

file_plate14 <- paste0(
  "/scratch/plapha/CP_240716 Batch effect test/output-newplatename/",
  "MyExpt_Cells.csv"
)


sce1 <- loadData(file_plate1)
sce2 <- loadData(file_plate2)
sce3 <- loadData(file_plate3)
sce4 <- loadData(file_plate4)
sce5 <- loadData(file_plate5)
sce6 <- loadData(file_plate6)
sce7 <- loadData(file_plate7)
sce8 <- loadData(file_plate8)
sce9 <- loadData(file_plate9)
sce10 <- loadData(file_plate10)
sce11 <- loadData(file_plate11)
sce12 <- loadData(file_plate12)
sce13 <- loadData(file_plate13)
sce14 <- loadData(file_plate14)

#sce <- cbind(sce1, sce2, sce3, sce4, sce5, sce6, sce7, sce8, sce9, sce10, sce11, sce12, sce13, sce14)

sce <- sce1
sce_111213 <- cbind(sce11, sce12, sce13)
sce_3456 <- cbind(sce3, sce4, sce5, sce6)
sce11 <- loadData(file_plate11)

```
Filter Sample overlapping between plate to test control -> validate permutation score method
#if score overlap on random curve = method is good

```{r}
# Get unique patient IDs from each experiment
patients_sce <- unique(sce$Patient)
patients_sce3 <- unique(sce3$Patient)
patients_sce4 <- unique(sce4$Patient)
patients_sce5 <- unique(sce5$Patient)
patients_sce6 <- unique(sce6$Patient)
patients_sce7 <- unique(sce7$Patient)
patients_sce8 <- unique(sce8$Patient)
patients_sce9 <- unique(sce9$Patient)
patients_sce10 <- unique(sce10$Patient)
patients_sce11 <- unique(sce11$Patients)
patients_sce12 <- unique(sce12$Patient)
patients_sce13 <- unique(sce13$Patient)
patients_sce14 <- unique(sce14$Patient)

patients_sce111213 <- unique(sce_111213$Patient)
patients_sce3456 <- unique(sce_3456$Patient)

```

Validate perturbation score with control cell from 2 different experiments, same Patient

```{r}
## 1. Find overlapping patients
common_patients <- intersect(patients_sce111213, patients_sce3456)

## 2. Filter both SCE objects for common patients
sce111213_filtered <- sce_111213[, sce_111213$Patient %in% common_patients]
sce3456_filtered <- sce_3456[, sce_3456$Patient %in% common_patients]

## 3. Add origin metadata to each object before merging
colData(sce111213_filtered)$source <- "sce_111213"
colData(sce3456_filtered)$source <- "sce_3456"

## 4. Merge the SCE objects
# Option 1: Simple column binding (if features are identical)
merged_sce <- cbind(sce111213_filtered, sce3456_filtered)

# filter to controls
cont_ids <- merged_sce$Treatment %in% c("control", "DMSO -IL1")
merged_sce <- merged_sce[, cont_ids]

# overwrite plate info - to tricl the code that they are from same plate
#merged_sce$Plate <- "fake"

```

Data transformation
#for validation

```{r}
# remove missing values in cells from merged sce
# don't need this step, the file was fakely merge as if they come from the same plate to be able to analyse

# define reference level
merged_sce$Treatment <- as.factor(merged_sce$Treatment)
merged_sce$Treatment <- relevel(merged_sce$Treatment, ref = "control")
 
# batch correct and normalize features
merged_sce <- normalizeExclude(merged_sce)
 
plotCellsPerImage(merged_sce)
merged_sce <- removeOutliers(merged_sce, min = 10, max = 300)


```

## Perturbation Score
#for validation

Compare perturbation score to null distribution.

```{r}
# overview of CellProfiler feature categories
#   include: Texture, Intensity, AreaShape, 
#   maybe include: RadialDistribution,  Neighbors, Location
#   exclude because many corner cases: Granularity, Correlation

sce_select <- merged_sce[str_detect(names(merged_sce), "Texture|Intensity|AreaShape"), ]

# perform PCA that will be used to calculate perturbation distances
sce_select <- runPCA(sce_select, exprs_values = "tfmfeatures")
plotReducedDim(sce_select, dimred = "PCA", colour_by = "Treatment")

# compute perturbation distances and relate to null distribution
treatmentEffect(sce_select, group = "Patient", treatment = "Treatment", 
                n_perm = 100, n_threads = 10)
treatmentEffect(sce_select, group = "Disease", treatment = "Treatment", 
                n_perm = 100, n_threads = 10)

```
#validation - Try split control within same plate into 'control' and 'test', check if it's has plate effect

```{r}
merged_sce <- sce111213_filtered
merged_sce$Treatment[merged_sce$Treatment == "control" , ]

## 1. Identify control cells
control_cells <- merged_sce[, merged_sce$Treatment == "control"]
cat("Original number of control cells:", ncol(control_cells), "\n")

## 2. Randomly split into two equal groups
set.seed(123) # For reproducibility
all_control_indices <- 1:ncol(control_cells)
split1_indices <- sample(all_control_indices, size = floor(ncol(control_cells)/2))
split2_indices <- setdiff(all_control_indices, split1_indices)

## 3. Create the two new SCE objects
control_split1 <- control_cells[, split1_indices]
control_split2 <- control_cells[, split2_indices]

## 4. Update Treatment names
control_split1$Treatment <- "control"  # Keep original name for group 1
control_split2$Treatment <- "test"     # New name for group 2

## 5. Verify the splits
cat("\nSplit results:\n")
cat("- Group 1 (control):", ncol(control_split1), "cells\n")
cat("- Group 2 (test):", ncol(control_split2), "cells\n")
cat("- Total:", ncol(control_split1) + ncol(control_split2), "cells (should match original)\n")

## 6. Create new combined SCE objects if needed
# Option A: New object with just these splits
new_sce_splits <- cbind(control_split1, control_split2)

# define reference level
new_sce_splits$Treatment <- as.factor(new_sce_splits$Treatment)
new_sce_splits$Treatment <- relevel(new_sce_splits$Treatment, ref = "control")
 
# batch correct and normalize features
new_sce_splits <- normalizeExclude(new_sce_splits)
 
plotCellsPerImage(new_sce_splits)
new_sce_splits <- removeOutliers(new_sce_splits, min = 20, max = 300)

##calculate pertubation score
sce_select <- new_sce_splits[str_detect(names(new_sce_splits), "Texture|Intensity|AreaShape"), ]

# perform PCA that will be used to calculate perturbation distances
sce_select <- runPCA(sce_select, exprs_values = "tfmfeatures")
plotReducedDim(sce_select, dimred = "PCA", colour_by = "Treatment")

# compute perturbation distances and relate to null distribution
treatmentEffect(sce_select, group = "Patient", treatment = "Treatment", 
                n_perm = 100, n_threads = 10)
treatmentEffect(sce_select, group = "Disease", treatment = "Treatment", 
                n_perm = 100, n_threads = 10)

```
Data transformation
#for actual experiment

```{r}
# remove missing values in cells
mat <- assay(sce, name = "features")
miss_cells <- apply(mat, 2, function(x) sum(is.na(x)))
cell_ids <- which(miss_cells > 0)
sce <- sce[, -cell_ids]
 
# define reference level
sce$Treatment <- as.factor(sce$Treatment)
sce$Treatment <- relevel(sce$Treatment, ref = "DMSO -IL1")
 
# batch correct and normalize features
sce <- normalizeExclude(sce)
 
plotCellsPerImage(sce)
sce <- removeOutliers(sce, min = 10, max = 200)


```

## Perturbation Score
#for accual experiment

Compare perturbation score to null distribution.

```{r}
# overview of CellProfiler feature categories
#   include: Texture, Intensity, AreaShape, 
#   maybe include: RadialDistribution,  Neighbors, Location
#   exclude because many corner cases: Granularity, Correlation

sce_select <- sce[str_detect(names(sce), "Texture|Intensity|AreaShape"), ]

# perform PCA that will be used to calculate perturbation distances
sce_select <- runPCA(sce_select, exprs_values = "tfmfeatures")
plotReducedDim(sce_select, dimred = "PCA", colour_by = "Treatment")

# compute perturbation distances and relate to null distribution
treatmentEffect(sce_select, group = "Patient", treatment = "Treatment", 
                n_perm = 100, n_threads = 10)
treatmentEffect(sce_select, group = "Disease", treatment = "Treatment", 
                n_perm = 100, n_threads = 10)
```





