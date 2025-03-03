library(readr)
library(ggplot2)
library(MASS)
library(dplyr)
library(tidyr)
library(stringr)
library(visdat)
library(factoextra)
library(FactoMineR)
library(uwot)
library(yardstick)
library(gridExtra)

library(cellpaintr)
library(scater)
library(scran)
library(bluster)
library(mbkmeans)
library(scuttle) 
library(viridis)
library(ComplexHeatmap)
library(circlize)




# plate 1

file_plate1 <- paste0(
  "/scratch/plapha/241212 citrulline shoulder/",
  "MyExpt_Cells.csv"
  )

file_plate2 <- paste0(
  "/scratch/plapha/241220 citrulline hand reimage/output/",
  "MyExpt_Cells.csv"
)

# load each plate into a separate bioconductor object
#plate1 <- read_csv(file_plate1)
#plate2 <- read_csv(file_plate2)

sce_plate1 <- loadData(file_plate1)
sce_plate2 <- loadData(file_plate2)

        #cells <- rbind(plate1, plate2)
sce <- cbind(sce_plate1, sce_plate2)

        # # number of cells per well
        # table(sce$Metadata_Patient)
        # ggplot(sce, aes(Metadata_Patient)) + geom_bar() + coord_flip()
        # 
        # 
        # # Filter dataset to include only Metadata_Treatment = "10%" and "0%"
        # filtered_cells <- cells %>% filter(Metadata_Treatment %in% c("10%", "0%"))
        # 
        # 
        # #Violin plot of intensity across cell types
        # 
        # ggplot(filtered_cells, aes(x = Metadata_Patient, y = Intensity_MeanIntensity_CorrCit, fill = Metadata_Treatment)) +
        #   geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(0.9)) +
        #   labs(title = "Protein Intensity Distribution for 10% vs. 0% Treatment",
        #        x = "Patient",
        #        y = "Protein Intensity",
        #        fill = "Treatment") +  # Legend for treatment comparison
        #   theme_minimal() +
        #   theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Data Transformation

plotCellsPerImage(sce)
sce <- removeOutliers(sce, min = 1, max = 100)
sce <- removeMissingValues(sce)
sce <- removeLowVariance(sce)
sce <- transformLogScale(sce)


#PCA plot
#Use `scater` for exploratory data analysis.

set.seed(38)
sce <- runPCA(sce, exprs_values = "tfmfeatures")
plotReducedDim(sce, dimred = "PCA", colour_by = "Patient")
    ggcells(sce, aes(x = PCA.1, y = PCA.2, color = Joint)) +
      geom_point() +
      theme_minimal() +
      labs(title = "PCA: PC1 vs PC2", x = "PC1", y = "PC2")
# plotReducedDim(sce, dimred = "PCA", colour_by = "Treatment")

###plot PCA correlation to interested PC --> PC2
    
    plotPCACor(sce, filter_by = 2, top = 30)
    
    
#see which feature contribute to PCA (PCA2)
    reducedDim(sce, "PCA")[,2]
    
#Extract PCA loadings
    pca_loadings <- attr(reducedDim(sce, "PCA"), "rotation") 
    
#Extract loadings for PCA 2 -corresponds to features (rows)
    pca2_loadings <- pca_loadings[, 2]
    
#Rank features by absolute contribution to PCA 2
    ranked_features <- sort(abs(pca2_loadings), decreasing = TRUE)
    
#Get the top 50 features
    top_50_features <- names(ranked_features)[1:50]
   
    print(top_50_features)
    
#plot loading to see most weighted loading (the elbow curve)
    
    # Get the feature names corresponding to the ranked loadings
    feature_names <- names(ranked_features) 
    
    # Create a data frame for plotting
    plot_data <- data.frame(
      Rank = seq_along(ranked_features),
      Loading = ranked_features,
      Feature = feature_names  # Add feature names as a column
    )
    
    
    top_n <- 10  # Label only the top 10 features
    # Plot the elbow curve with feature names
    
    library(ggrepel)
    
    ggplot(plot_data, aes(x = Rank, y = Loading)) +
      geom_line() +
      geom_point() +
      geom_text_repel(
        data = subset(plot_data, Rank <= top_n),  # Label only top N features
        aes(label = Feature),
        hjust = -0.2, vjust = 0.5, size = 3, direction = "y"  # Adjust hjust and vjust as needed
      ) +
      labs(
        title = "Elbow Curve of Loading Values for PC2",
        x = "Rank of Features",
        y = "Absolute Loading Value"
      ) +
      theme_minimal()
    
    
#tSNE
set.seed(23)
sce <- runTSNE(sce, perplexity = 80, exprs_values = "tfmfeatures",
               num_threads = 16) # set to the number of cores for speed up
plotReducedDim(sce, dimred = "TSNE", colour_by = "Patient")
# plotReducedDim(sce, dimred = "TSNE", colour_by = "Treatment")
# plotReducedDim(sce, dimred = "TSNE", colour_by = "Joint")

#UMAP
set.seed(23)
sce <- runUMAP(sce, exprs_values = "tfmfeatures")
plotReducedDim(sce, dimred = "UMAP", colour_by = "Patient")

# plotReducedDim(sce, dimred = "UMAP", colour_by = "Treatment")
# plotReducedDim(sce, dimred = "UMAP", colour_by = "Joint")
# plotReducedDim(sce, dimred = "UMAP", colour_by = "Disease")

#Custom plots using `ggcells` function.
# ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2, colour = Treatment)) +
#   geom_point(alpha = 0.5) +
#   facet_wrap(~Patient) +
#   theme_minimal()

##Plot each patient facet on UMAP

  #Extract UMAP Coordinates and Metadata
umap_coords <- reducedDim(sce, "UMAP")
metadata <- colData(sce)

  #Combine UMAP Coordinates and Metadata
umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Patient = metadata$Patient
)

  #Plot UMAP with Facets
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Patient), size = 1, alpha = 0.7) +
  facet_wrap(~ Patient) +  # Facet by Patient
  theme_minimal() +
  labs(title = "UMAP Colored by Patient", x = "UMAP1", y = "UMAP2")


###############################################################

#Clustering from Christof

nn.clusters <- clusterCells(sce, use.dimred = "PCA", 
                            BLUSPARAM = SNNGraphParam(k = 10, type = "rank", cluster.fun = "walktrap"))
table(nn.clusters)
colLabels(sce) <- nn.clusters
plotReducedDim(sce, "UMAP", colour_by = "label")

#Faster clustering: 
      ##before was center = 1000, k=10
      #reduced to 2 clusters
set.seed(23) 
kgraph.clusters <- clusterCells(sce, use.dimred = "PCA",
                                BLUSPARAM = TwoStepParam(
                                  first = KmeansParam(centers = 1000),
                                  second = SNNGraphParam(k = 10, type = "rank", cluster.fun = "walktrap")
                                )
)

table(kgraph.clusters)
colLabels(sce) <- kgraph.clusters
plotReducedDim(sce, "PCA", colour_by = "label")
plotReducedDim(sce, "UMAP", colour_by = "label")


#####this one work########use Louvain for clustering ##only seed from 36-39
set.seed(38)
kgraph.clusters <- clusterCells(sce, use.dimred = "PCA",
                                BLUSPARAM = TwoStepParam(
                                  first = KmeansParam(centers = 10),
                                  second = SNNGraphParam(k = 4, type = "rank", cluster.fun = "louvain") # Try Louvain instead
                                )
)


table(kgraph.clusters)
colLabels(sce) <- kgraph.clusters
plotReducedDim(sce, "PCA", colour_by = "label")
plotReducedDim(sce, "UMAP", colour_by = "label")


# set UMAP coordinate as a ref for clustering,PLUS set threshold on UMAP1 at -3
umap_coords <- reducedDim(sce, "UMAP")
km.clusters <- ifelse(umap_coords[,1] < -3, "Cluster 2", "Cluster 1")
km.clusters <- factor(km.clusters, levels = c("Cluster 1", "Cluster 2"))
colLabels(sce) <- km.clusters
plotReducedDim(sce, "UMAP", colour_by = "label")

#################################################################

#Marker gene detection by pair-wise comparison

#make sure 'logcounts'column is in sce, ready for scoreMarkers function

assay(sce, "logcounts") <- assay(sce, "tfmfeatures")
marker.info <- scoreMarkers(sce, colLabels(sce))
marker.info

colnames(marker.info[["1"]]) # statistics for cluster 1.


#see candidate markers from cluster

chosen <- marker.info[["2"]] #select cluster
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])

plotExpression(sce, features=head(rownames(ordered)), 
               x="label", colour_by="label")

#Cohen's d (standardized log-fold change) show effect size of each feature in cluster
cohen.only <- chosen[,grepl("logFC.cohen", colnames(chosen))]
cohen.only[order(cohen.only$mean.logFC.cohen,decreasing=TRUE),]


#summarized pairwise effect --> use 'median' with Cohen from previous step
chosen <- marker.info[["1"]] # using another cluster, for some variety.
ordered <- chosen[order(chosen$median.logFC.cohen,decreasing=TRUE),]
head(ordered[,1:4]) 

plotExpression(sce, features=head(rownames(ordered)), 
               x="label", colour_by="label")

#summarized by ranking
ordered <- chosen[order(chosen$rank.logFC.cohen),]
top.ranked <- ordered[ordered$rank.logFC.cohen <= 5,]
rownames(top.ranked)

plotGroupedHeatmap(sce, features=rownames(top.ranked), group="label", 
                   center=TRUE, zlim=c(-3, 3))

# Omitting the decreasing=TRUE to focus on negative effects(down regulated features)
ordered <- chosen[order(chosen$median.logFC.cohen),1:4]
head(ordered)

###Using a log-fold change threshold

#cluster1
marker.info.lfc <- scoreMarkers(sce, colLabels(sce), lfc=1)
chosen1 <- marker.info.lfc[["1"]] # another cluster for some variety.
chosen1 <- chosen2[order(chosen2$mean.AUC, decreasing=TRUE),]
chosen1[,c("self.average", "other.average", "mean.AUC")]


plotDots(sce, rownames(chosen1)[1:20], group="label")

#cluster2
marker.info.lfc <- scoreMarkers(sce, colLabels(sce), lfc=1)
chosen2 <- marker.info.lfc[["2"]] # another cluster for some variety.
chosen2 <- chosen2[order(chosen2$mean.AUC, decreasing=TRUE),]
chosen2[,c("self.average", "other.average", "mean.AUC")]


plotDots(sce, rownames(chosen2)[1:20], group="label")


##########plot heatmap from differential features

library(viridis)
library(ComplexHeatmap)

#Step1 Extract metadata for annotation
metadata <- colData(sce)[, c("Gender", "Joint", "Patient", "Treatment", "Disease")]  # Replace with your metadata columns

#Step2 Automatically generate colors for each category
gender_colors <- viridis(length(unique(metadata$Gender)))
names(gender_colors) <- unique(metadata$Gender)

joint_colors <- viridis(length(unique(metadata$Joint)))
names(joint_colors) <- unique(metadata$Joint)

patient_colors <- viridis(length(unique(metadata$Patient)))
names(patient_colors) <- unique(metadata$Patient)

treatment_colors <- viridis(length(unique(metadata$Treatment)))
names(treatment_colors) <- unique(metadata$Treatment)

disease_colors <- viridis(length(unique(metadata$Disease)))
names(disease_colors) <- unique(metadata$Disease) 

#Step3 Create HeatmapAnnotation
ha <- HeatmapAnnotation(
  df = metadata,
  col = list(
    Joint = joint_colors,
    Gender = gender_colors,
    Patient = patient_colors,
    Treatment = treatment_colors,
    Disease = disease_colors
  )
)

# Step 4: Extract differentially expressed features based on LFC and mean.AUC
marker.info.lfc <- scoreMarkers(sce, colLabels(sce), lfc = 1)

# Initialize an empty list to store selected DEGs
selected_genes <- list()

# Loop through clusters 1 to 7
for (cluster in 1:7) {
  # Extract DEGs for the current cluster
  chosen <- marker.info.lfc[[as.character(cluster)]]
  
  # Order DEGs by mean.AUC in descending order
  chosen <- chosen[order(chosen$mean.AUC, decreasing = TRUE), ]
  
  # Store the top DEGs in the list
  selected_genes[[paste0("Cluster_", cluster)]] <- chosen
}

# Combine top DEGs from all clusters, select only 50 features (Ensure top_genes is a character vector)
      # #top_genes <- unique(unlist(lapply(selected_genes, function(x) {
      #  # if (length(x) > 0) return(x[1:min(50, length(x))]) else return(NULL)
      # })))
top_genes <- unique(unlist(lapply(selected_genes, function(x) rownames(x)[1:50])))

valid_genes <- top_genes[top_genes %in% rownames(logcounts_matrix)]

logcounts_subset_top <- logcounts_matrix[valid_genes, ]

# Print subset size
cat("Final subset size:", dim(logcounts_subset_top), "\n")  




# Step 5: Extract logcounts assay and subset for top DEGs
logcounts_matrix <- logcounts(sce)  # Extract logcounts assay
logcounts_subset_top <- logcounts_matrix[top_genes, ]  # Subset for top DEGs

# Step 6: Plot the heatmap
Heatmap(
  logcounts_subset_top,
  name = "Expression",
  top_annotation = ha,
  show_row_names = FALSE,  # Show row names (gene names)
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = metadata$Cluster  
)


#######################Heatmap from already cluster data

#make sure 'logcounts'column is in sce, ready for scoreMarkers function

assay(sce, "logcounts") <- assay(sce, "tfmfeatures")

# Step 1: Identify Top DEGs per Cluster
marker.info <- scoreMarkers(sce, groups = colLabels(sce))  # Identify marker genes
selected_genes <- lapply(marker.info, function(x) rownames(x[order(x$mean.AUC, decreasing = TRUE), ])[1:50])  # Top 50 genes per cluster
top_genes <- unique(unlist(selected_genes))  # Get unique list of top genes


# Step 2: Subset Logcounts Matrix for Selected Genes
logcounts_matrix <- logcounts(sce)  # Extract log-normalized counts
valid_genes <- intersect(rownames(logcounts_matrix), top_genes)  # Keep only matching genes
logcounts_subset_top <- logcounts_matrix[valid_genes, ]  # Subset for heatmap

#####plot in order of CP features
    # Step 2.1: Extract Feature Names
    feature_names <- rownames(rowData(sce))
    
    ## Step 2.2: Define Feature Categories
    categories <- list(
      Intensity = grep("Intensity", feature_names, value = TRUE, ignore.case = TRUE),
      Texture   = grep("Texture", feature_names, value = TRUE, ignore.case = TRUE),
      Shape     = grep("Shape|FormFactor|Eccentricity|Perimeter", feature_names, value = TRUE, ignore.case = TRUE),
      Area      = grep("Area", feature_names, value = TRUE, ignore.case = TRUE)
    )

    # Step 2.3: Identify Features Not in Main Categories
    all_categorized_features <- unlist(categories)
    categories$Other <- setdiff(feature_names, all_categorized_features)
    
    
# Step 3: Create Metadata for Annotation
metadata <- colData(sce)[, c("Joint", "Patient", "Treatment", "Disease")]


# Step 4: Generate Annotation Colors
joint_colors <- viridis(length(unique(metadata$Joint)))
names(joint_colors) <- unique(metadata$Joint)

patient_colors <- viridis(length(unique(metadata$Patient)))
names(patient_colors) <- unique(metadata$Patient)

treatment_colors <- viridis(length(unique(metadata$Treatment)))
names(treatment_colors) <- unique(metadata$Treatment)

disease_colors <- viridis(length(unique(metadata$Disease)))
names(disease_colors) <- unique(metadata$Disease)

# Step 5: Create Heatmap Annotation
ha <- HeatmapAnnotation(
  df = metadata,
  col = list(
    Joint = joint_colors,
    Patient = patient_colors,
    Treatment = treatment_colors,
    Disease = disease_colors
  )
)

# Step 6: Plot the Heatmap
Heatmap(
  logcounts_subset_top,
  name = "Expression",
  top_annotation = ha,
  show_row_names = FALSE,  # Show gene names
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = colLabels(sce)  # Split by cluster
)

        # Step 8: Plot the Heatmap based on CP features categories
        
    categories <- lapply(categories, unique)
    row_split <- factor(valid_genes, levels = unique(unlist(categories)))

          Heatmap(
            logcounts_subset_top,
            name = "Expression",
            top_annotation = ha,
            show_row_names = FALSE,  # Hide row names for readability
            show_column_names = FALSE,
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            column_split = colLabels(sce),  # Split by cluster
            row_split = row_split  # Use the corrected factor for row splitting
          )

###################################Check and Group Features into categories#######
 #make sure 'logcounts'column is in sce, ready for scoreMarkers function
          
 assay(sce, "logcounts") <- assay(sce, "tfmfeatures")
          
#Step 1: Extract Feature Names
feature_names <- rownames(rowData(sce))

#Step 2: Define Feature Categories Based on Keywords
categories <- list(
  Intensity = grep("Intensity", feature_names, value = TRUE, ignore.case = TRUE),
  Texture   = grep("Texture", feature_names, value = TRUE, ignore.case = TRUE),
  Shape     = grep("Shape|FormFactor|Eccentricity|Perimeter", feature_names, value = TRUE, ignore.case = TRUE),
  Area      = grep("Area", feature_names, value = TRUE, ignore.case = TRUE)
)

# Step 3: Identify Features Not in the Main Categories
all_categorized_features <- unlist(categories)
categories$Other <- setdiff(feature_names, all_categorized_features)

# # Step 4: Print Summary of Features in Each Category
# for (cat in names(categories)) {
#   cat("\n🔹", cat, "Features:", length(categories[[cat]]), "\n")
#   print(head(categories[[cat]], 10))  # Show first 10 features in each category
# }

##### now plot heatmap from the categories

# Step 1: Identify Top DEGs per Cluster
marker.info <- scoreMarkers(sce, groups = colLabels(sce))  # Identify marker genes
selected_genes <- lapply(marker.info, function(x) rownames(x[order(x$mean.AUC, decreasing = TRUE), ])[1:50])  # Top 50 genes per cluster
top_genes <- unique(unlist(selected_genes)) 

# Step 2: Subset Logcounts Matrix for Selected Genes
logcounts_matrix <- logcounts(sce)  # Extract log-normalized counts
valid_genes <- intersect(rownames(logcounts_matrix), top_genes)  # Keep only matching genes
logcounts_subset_top <- logcounts_matrix[valid_genes, ] 

# Step 2.1: Extract Feature Names
feature_names <- rownames(rowData(sce))

# Step 2.2: Define Feature Categories
categories <- list(
  Intensity = grep("Intensity", feature_names, value = TRUE, ignore.case = TRUE),
  Texture   = grep("Texture", feature_names, value = TRUE, ignore.case = TRUE),
  Shape     = grep("Shape|FormFactor|Eccentricity|Perimeter", feature_names, value = TRUE, ignore.case = TRUE),
  Area      = grep("Area", feature_names, value = TRUE, ignore.case = TRUE)
)

# Step 2.3: Identify Features Not in Main Categories
all_categorized_features <- unlist(categories)
categories$Other <- setdiff(feature_names, all_categorized_features)

# Step 3: Create Metadata for Annotation
metadata <- colData(sce)[, c("Joint", "Patient", "Treatment", "Disease")]

# Step 4: Generate Annotation Colors
joint_colors <- viridis(length(unique(metadata$Joint)))
names(joint_colors) <- unique(metadata$Joint)

patient_colors <- viridis(length(unique(metadata$Patient)))
names(patient_colors) <- unique(metadata$Patient)

treatment_colors <- viridis(length(unique(metadata$Treatment)))
names(treatment_colors) <- unique(metadata$Treatment)

disease_colors <- viridis(length(unique(metadata$Disease)))
names(disease_colors) <- unique(metadata$Disease)

# Step 5: Create Heatmap Annotation
ha <- HeatmapAnnotation(
  df = metadata,
  col = list(
    Joint = joint_colors,
    Patient = patient_colors,
    Treatment = treatment_colors,
    Disease = disease_colors
  )
)

# Step 6: Plot Heatmaps for Each Category#### save files in working directory for R to run multiple heatmap
for (category_name in names(categories)) {
  category_genes <- categories[[category_name]]
  valid_category_genes <- intersect(valid_genes, category_genes)
  
  if (length(valid_category_genes) > 0) {
    logcounts_subset_category <- logcounts_subset_top[valid_category_genes, ]
    
    # Save heatmap to a file
    png(paste0("heatmap_", category_name, ".png"), width = 800, height = 600)
    draw(
      Heatmap(
        logcounts_subset_category,
        name = "Expression",
        top_annotation = ha,
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_split = colLabels(sce),
        row_split = factor(valid_category_genes, levels = valid_category_genes)
      )
    )
    dev.off()  # Close the file device
  }
}

#####only plot sub-category

##intensity
intensity_genes <- categories$Intensity
valid_intensity_genes <- intersect(valid_genes, intensity_genes)
logcounts_subset_intensity <- logcounts_subset_top[valid_intensity_genes, ]


# Define a color scale

library(circlize)


col_fun <- colorRamp2(
  c(min(logcounts_subset_intensity), max(logcounts_subset_intensity)),  # Breaks (two values)
  viridis(2)  # Only two colors (start and end)
    )

#add cluster number to heatmap
cluster_numbers <- km.clusters
cluster_colors <- c("Cluster 1" = "steelblue2", "Cluster 2" = "darkorange1") 

# Column annotations
ha <- HeatmapAnnotation(
  df = metadata,
  col = list(
    Joint = joint_colors,
    Patient = patient_colors,
    Treatment = treatment_colors,
    Disease = disease_colors,
    Cluster = cluster_colors  # Add cluster colors
  ),
  Cluster = cluster_numbers  # Add cluster numbers as an annotation
)



# Generate heatmap
png("intensity_heatmap_2cluster_noLocation.png", width = 1200, height = 800)
draw(
  Heatmap(
    logcounts_subset_intensity,
    name = "Expression",
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_split = colLabels(sce),
    col = col_fun,
    column_title = "Intensity Sub-Category Heatmap",
    width = unit(10, "cm"),
    height = unit(15, "cm")
  )
)
dev.off()


##Shape
shape_genes <- categories$Shape
valid_shape_genes <- intersect(valid_genes, shape_genes)
logcounts_subset_shape <- logcounts_subset_top[valid_shape_genes, ]


# Define a color scale

library(circlize)


col_fun <- colorRamp2(
  c(min(logcounts_subset_shape), max(logcounts_subset_shape)),  # Breaks (two values)
  viridis(2)  # Only two colors (start and end)
)

#add cluster number to heatmap
cluster_numbers <- km.clusters
cluster_colors <- c("Cluster 1" = "steelblue2", "Cluster 2" = "darkorange1") 

# Column annotations
ha <- HeatmapAnnotation(
  df = metadata,
  col = list(
    Joint = joint_colors,
    Patient = patient_colors,
    Treatment = treatment_colors,
    Disease = disease_colors,
    Cluster = cluster_colors  # Add cluster colors
  ),
  Cluster = cluster_numbers  # Add cluster numbers as an annotation
)

# Column annotations
ha <- HeatmapAnnotation(
  df = metadata,
  col = list(
    Joint = joint_colors,
    Patient = patient_colors,
    Treatment = treatment_colors,
    Disease = disease_colors
  )
)


# Generate heatmap
png("shape_heatmap_2cluster.png", width = 1200, height = 800)
draw(
  Heatmap(
    logcounts_subset_shape,
    name = "Expression",
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_split = colLabels(sce),
    col = col_fun,
    column_title = "Shape Sub-Category Heatmap",
    width = unit(10, "cm"),
    height = unit(15, "cm")
  )
)
dev.off()


##other features
other_genes <- categories$Other
valid_other_genes <- intersect(valid_genes, other_genes)
logcounts_subset_other <- logcounts_subset_top[valid_other_genes, ]


# Define a color scale

library(circlize)


col_fun <- colorRamp2(
  c(min(logcounts_subset_other), max(logcounts_subset_other)),  # Breaks (two values)
  viridis(2)  # Only two colors (start and end)
)


# Column annotations
ha <- HeatmapAnnotation(
  df = metadata,
  col = list(
    Joint = joint_colors,
    Patient = patient_colors,
    Treatment = treatment_colors,
    Disease = disease_colors
  )
)


# Generate heatmap
png("other_heatmap_2cluster.png", width = 1200, height = 800)
draw(
  Heatmap(
    logcounts_subset_other,
    name = "Expression",
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_split = colLabels(sce),
    col = col_fun,
    column_title = "other Sub-Category Heatmap",
    width = unit(10, "cm"),
    height = unit(15, "cm")
  )
)
dev.off()


#######Area
area_genes <- categories$Area
valid_area_genes <- intersect(valid_genes, area_genes)
logcounts_subset_area <- logcounts_subset_top[valid_area_genes, ]


# Define a color scale

library(circlize)


col_fun <- colorRamp2(
  c(min(logcounts_subset_area), max(logcounts_subset_area)),  # Breaks (two values)
  viridis(2)  # Only two colors (start and end)
)


# Column annotations
ha <- HeatmapAnnotation(
  df = metadata,
  col = list(
    Joint = joint_colors,
    Patient = patient_colors,
    Treatment = treatment_colors,
    Disease = disease_colors
  )
)


# Generate heatmap
png("area_heatmap_2cluster.png", width = 1200, height = 800)
draw(
  Heatmap(
    logcounts_subset_area,
    name = "Expression",
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_split = colLabels(sce),
    col = col_fun,
    column_title = "Area Sub-Category Heatmap",
    width = unit(10, "cm"),
    height = unit(15, "cm")
  )
)
dev.off()

###############
####################
# Filter Location feature out

# Identify features that start with "Location_"
location_features <- grep("^Location_", rownames(sce), value = TRUE)

# Remove these features from the sce object
sce <- sce[!rownames(sce) %in% location_features, ]

#PCA
set.seed(38)
sce <- runPCA(sce, exprs_values = "tfmfeatures")
plotReducedDim(sce, dimred = "PCA", colour_by = "Joint")

#UMAP
set.seed(38)
sce <- runUMAP(sce, exprs_values = "tfmfeatures")
plotReducedDim(sce, dimred = "UMAP", colour_by = "Patient")

# set UMAP coordinate as a ref for clustering,PLUS set threshold on UMAP1 at 5
umap_coords <- reducedDim(sce, "UMAP")
km.clusters <- ifelse(umap_coords[,1] < 5, "Cluster 2", "Cluster 1")
km.clusters <- factor(km.clusters, levels = c("Cluster 1", "Cluster 2"))
colLabels(sce) <- km.clusters
plotReducedDim(sce, "PCA", colour_by = "label")

metadata <- as.data.frame(colData(sce))
cluster_labels <- colLabels(sce)

ggplot(metadata, aes(x = factor(cluster_labels), fill = Joint)) +
  geom_bar(position = "fill") +  # Normalize each bar to 100%
  scale_y_continuous(labels = scales::percent_format()) +  # Convert y-axis to percentages
  labs(title = "Joint Distribution by Cluster", 
       x = "Cluster", 
       y = "Percentage") +
  theme_minimal()


#Check if specific metadata categories are associated with differentially expressed features.
feature_of_interest <- logcounts(sce)["Intensity_IntegratedIntensity_CorrCit", ]  # Replace with your feature
metadata$feature_of_interest <- feature_of_interest

# Logistic regression to test association
logistic_model <- glm(feature_of_interest ~ Disease, data = metadata, family = binomial)
summary(logistic_model)

linear_model <- lm(feature_of_interest ~ Disease, data = metadata)
summary(linear_model)


##To separate the feature data into different channels (e.g., CorrCit, CorrCyto, and CorrBlue) 
#and run UMAP for each of them

# Extract feature names
feature_names <- rownames(rowData(sce))

# Subset features for each channel
corr_cit_features <- feature_names[grep("CorrCit", feature_names)]
corr_cyto_features <- feature_names[grep("CorrCyto", feature_names)]
corr_blue_features <- feature_names[grep("CorrBlue", feature_names)]

# Extract data for each channel

assay(sce, "logcounts") <- assay(sce, "tfmfeatures")

corr_cit_data <- logcounts(sce)[corr_cit_features, ]
corr_cyto_data <- logcounts(sce)[corr_cyto_features, ]
corr_blue_data <- logcounts(sce)[corr_blue_features, ]

#Create New SingleCellExperiment Objects for Each Channel
library(SingleCellExperiment)

# Create SingleCellExperiment objects for each channel
sce_corr_cit <- SingleCellExperiment(assays = list(logcounts = corr_cit_data), colData = colData(sce))
sce_corr_cyto <- SingleCellExperiment(assays = list(logcounts = corr_cyto_data), colData = colData(sce))
sce_corr_blue <- SingleCellExperiment(assays = list(logcounts = corr_blue_data), colData = colData(sce))

#Run UMAP on each SingleCellExperiment object
set.seed(38)

# Run UMAP for CorrCit
sce_corr_cit <- runUMAP(sce_corr_cit, exprs_values = "logcounts")
plotReducedDim(sce_corr_cit, dimred = "UMAP", colour_by = "Treatment") +
  ggtitle("UMAP for CorrCit Features")

# umap_coords <- reducedDim(sce_corr_cit, "UMAP")
# plot_data <- data.frame(
#   UMAP1 = umap_coords[, 1],
#   UMAP2 = umap_coords[, 2],
#   Treatment = colData(sce_corr_cit)$Treatment,  # Ensure "Treatment" exists in colData
#   Patient = colData(sce_corr_cit)$Patient      # Ensure "Patient" exists in colData
# )
# 
# ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Treatment)) +
#   geom_point(size = 1, alpha = 0.2) +
#   theme_minimal() +
#   labs(
#     title = "UMAP for CorrCit Features",
#     x = "UMAP1",
#     y = "UMAP2"
#   ) +
#   facet_wrap(~Patient) +  # Facet by Patient
#   theme(legend.position = "top")


# Run UMAP for CorrCyto
sce_corr_cyto <- runUMAP(sce_corr_cyto, exprs_values = "logcounts")
plotReducedDim(sce_corr_cyto, dimred = "UMAP", colour_by = "Treatment") +
  ggtitle("UMAP for CorrCyto Features")


# Run UMAP for CorrBlue
sce_corr_blue <- runUMAP(sce_corr_blue, exprs_values = "logcounts")
plotReducedDim(sce_corr_blue, dimred = "UMAP", colour_by = "Treatment") +
  ggtitle("UMAP for CorrBlue Features")

#Run PCA for each channel
sce_corr_cit <- runPCA(sce_corr_cit, exprs_values = "logcounts")
plotReducedDim(sce_corr_cit, dimred = "PCA", colour_by = "Treatment") + ggtitle("PCA for CorrCit Features")

sce_corr_cyto <- runPCA(sce_corr_cyto, exprs_values = "logcounts")
plotReducedDim(sce_corr_cyto, dimred = "PCA", colour_by = "Treatment") + ggtitle("PCA for CorrCyto Features")

sce_corr_blue <- runPCA(sce_corr_blue, exprs_values = "logcounts")
plotReducedDim(sce_corr_blue, dimred = "PCA", colour_by = "Treatment") + ggtitle("PCA for CorrBlue Features")


##Plot integrated intensity on UMAP to see if treatment arrest cells

#Extract the UMAP embeddings for CorrBlue
umap_coords <- reducedDim(sce_corr_blue, "UMAP")

#Extract the feature values for Intensity_Integratedintensity_CorrBlue
intensity_feature <- logcounts(sce_corr_blue)["Intensity_IntegratedIntensity_CorrBlue", ]

#Create a data frame for plotting
treatment_info <- colData(sce_corr_blue)$Treatment
plot_data <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Intensity = intensity_feature,
  Treatment = treatment_info 
)

#Plot the feature overlay on UMAP
ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Intensity)) +
  geom_point(size = 1, alpha = 0.2) +
  scale_color_viridis_c(name = "Intensity") +  # Use a color scale (e.g., viridis)
  theme_minimal() +
  labs(
    title = "UMAP for CorrBlue with Intensity_IntegratedIntensity_CorrBlue Overlay",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  facet_wrap(~Treatment)

#Plot histogram of Intensity_IntegratedIntensity_CorrBlue 

treatment_info <- colData(sce_corr_blue)$Treatment
patient_info <- colData(sce_corr_blue)$Patient

plot_data <- data.frame(
  Intensity = intensity_feature,
  Treatment = treatment_info,
  Patient = patient_info
)

ggplot(plot_data, aes(x = Intensity)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Histogram of Intensity_Integratedintensity_CorrBlue",
    x = "Intensity",
    y = "Count"
  ) +
  facet_wrap(~Treatment) 


#plot line histogram instead of bar and overlay on-top of each other

ggplot(plot_data, aes(x = Intensity, color = Treatment, fill = Treatment)) +
  geom_density(alpha = 0.3, adjust = 1.2) +  # Density plot with transparency
  theme_minimal() +
  labs(
    title = "Density Plot of Intensity_Integratedintensity_CorrBlue",
    x = "Intensity",
    y = "Density"
  ) +
  scale_fill_manual(values = c("red", "blue", "green", "purple")) +  # Custom colors for treatments
  scale_color_manual(values = c("red", "blue", "green", "purple")) +  
  theme(legend.position = "top") +   # Move legend to top
  facet_wrap(~Patient) 

#Clustering for sce_corr_cit 
#Faster clustering: 

set.seed(36) 
kgraph.clusters <- clusterCells(sce_corr_cit, use.dimred = "UMAP",
                                BLUSPARAM = TwoStepParam(
                                  first = KmeansParam(centers = 100),
                                  second = SNNGraphParam(k = 5, type = "rank", cluster.fun = "walktrap")
                                )
)

table(kgraph.clusters)
colLabels(sce_corr_cit) <- kgraph.clusters

plotReducedDim(sce_corr_cit, "UMAP", colour_by = "label")

#plot heatmap from cluster
# Step 1: Identify Top DEGs per Cluster
marker.info <- scoreMarkers(sce_corr_cit, groups = colLabels(sce_corr_cit))  # Identify marker genes
selected_genes <- lapply(marker.info, function(x) rownames(x[order(x$mean.AUC, decreasing = TRUE), ])[1:50])  # Top 50 genes per cluster
top_genes <- unique(unlist(selected_genes))  # Get unique list of top genes


# Step 2: Subset Logcounts Matrix for Selected Genes
logcounts_matrix <- logcounts(sce_corr_cit)  # Extract log-normalized counts
valid_genes <- intersect(rownames(logcounts_matrix), top_genes)  # Keep only matching genes
logcounts_subset_top <- logcounts_matrix[valid_genes, ] 

# Step 3: Create Metadata for Annotation
metadata <- colData(sce_corr_cit)[, c("Joint", "Patient", "Treatment", "Disease")]


# Step 4: Generate Annotation Colors
joint_colors <- viridis(length(unique(metadata$Joint)))
names(joint_colors) <- unique(metadata$Joint)

patient_colors <- viridis(length(unique(metadata$Patient)))
names(patient_colors) <- unique(metadata$Patient)

treatment_colors <- viridis(length(unique(metadata$Treatment)))
names(treatment_colors) <- unique(metadata$Treatment)

disease_colors <- viridis(length(unique(metadata$Disease)))
names(disease_colors) <- unique(metadata$Disease)

# Step 5: Create Heatmap Annotation
ha <- HeatmapAnnotation(
  df = metadata,
  col = list(
    Joint = joint_colors,
    Patient = patient_colors,
    Treatment = treatment_colors,
    Disease = disease_colors
  )
)

# Step 6: Plot the Heatmap
Heatmap(
  logcounts_subset_top,
  name = "Expression",
  top_annotation = ha,
  show_row_names = TRUE,  # Show gene names
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = colLabels(sce_corr_cit)  # Split by cluster
)




########Load new csv file from Nucleaus

file_plate1 <- paste0(
  "/scratch/plapha/241212 citrulline shoulder/",
  "MyExpt_Nuclei.csv"
)

file_plate2 <- paste0(
  "/scratch/plapha/241220 citrulline hand reimage/output/",
  "MyExpt_Nuclei.csv"
)

sce_plate1 <- loadData(file_plate1)
sce_plate2 <- loadData(file_plate2)

sce <- cbind(sce_plate1, sce_plate2)

######################################

###From PCA -- want to see which feature contribute to PC that has separation
## in this data set - it is PC2

#Extract loadings for PCA 2 -corresponds to features (rows)
pca2_loadings <- pca_loadings[, 2]

# Step 1: Create the top 50 features table
pca2_loadings_df <- data.frame(
  Feature = rownames(pca_loadings),
  Loading = pca2_loadings
)

# Sort by absolute loadings
pca2_loadings_df <- pca2_loadings_df[order(abs(pca2_loadings_df$Loading), decreasing = TRUE), ]

# Get the top 50 features
top_50_features_df <- pca2_loadings_df[1:50, ]

# Step 2: Export the table to a CSV file
write.csv(top_50_features_df, file = "top_50_features_PCA2.csv", row.names = FALSE)



###### now track back to which Metadata drive PC2

#Correlate Metadata with PCA2 Scores - For Categorical Metadata (e.g., Patient, Treatment, Joint):
#Use ANOVA or Kruskal-Wallis 

#If the sample size is too large (e.g., > 5000), use a Q-Q plot to visually assess normality

# Q-Q plot for PCA2 scores
qqnorm(pca2_scores)
qqline(pca2_scores, col = "red")

#If the points fall approximately along the red line, the data is normally distributed. - use ANOVA

#If the points deviate significantly from the line, the data is not normally distributed. - use Kruskal-Wallis 


#Extract the PCA2 scores -correspond to samples (column)
pca2_scores <- reducedDim(sce, "PCA")[, 2]

# Step 2: Combine PCA2 scores with metadata
metadata <- colData(sce)
pca2_metadata_df <- data.frame(
  PC2 = pca2_scores,
  Joint = metadata$Joint,
  Treatment = metadata$Treatment,
  Patient = metadata$Patient,
  Disease = metadata$Disease
)

# Step 3: Run Kruskal-Wallis test
kruskal_result <- kruskal.test(PC2 ~ Patient, data = pca2_metadata_df)
print(kruskal_result)


#To visualize how PCA2 scores vary across Joint
# Boxplot for Joint
ggplot(pca2_metadata_df, aes(x = Disease, y = PC2, fill = Disease)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "PCA2 Scores by Disease",
    x = "Disease",
    y = "PCA2 Scores"
  )+
  theme(
    plot.title = element_text(size = 20),  # Increase title font size
    axis.title.x = element_text(size = 16),  # Increase x-axis label font size
    axis.title.y = element_text(size = 16),  # Increase y-axis label font size
    axis.text.x = element_text(size = 14),  # Increase x-axis tick labels font size
    axis.text.y = element_text(size = 14),  # Increase y-axis tick labels font size
    legend.title = element_text(size = 16),  # Increase legend title font size
    legend.text = element_text(size = 14)  # Increase legend text font size
  )
