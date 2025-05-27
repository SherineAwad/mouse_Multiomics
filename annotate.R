library(ArchR)
library(dplyr)
library(ggplot2)
library(viridis)

# Setup
addArchRThreads(threads = 4)
addArchRGenome("mm10")
setwd("/nfs/turbo/umms-thahoang/sherine/mouse_Multiomics")


args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]
annotation_file <- args[2] 
proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)


annotation_df <- read.csv(annotation_file, stringsAsFactors = FALSE)  # or read.delim/read_tsv as needed

# Convert to named vector
cluster_to_celltype <- setNames(annotation_df$CellType, annotation_df$Cluster)

# Add to ArchR project metadata
proj_ALL$CellType <- cluster_to_celltype[proj_ALL$Clusters_Combined]


###
# ================================
# Plot and Save UMAP with Cell Type Labels
# ================================

figure_name <- paste0(project_name, "_CellTypeUMAP_annotated.png")

p <- plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "CellType",
  embedding = "UMAP_Combined",
  labelAsFactors = FALSE,
  labelMeans = TRUE  # This adds the cell type labels at cluster centers
)

ggsave(filename = figure_name, plot = p, width = 6, height = 5, dpi = 300)

# ================================
# Save Annotated Project
# ================================

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = paste0(project_name, "_annotated"), load = FALSE)

