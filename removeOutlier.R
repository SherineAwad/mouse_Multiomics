library(ArchR)
library(dplyr)
library(ggplot2)
library(viridis)

# Setup
addArchRThreads(threads = 4)
addArchRGenome("mm10")
setwd("/nfs/turbo/umms-thahoang/sherine/mouse_Multiomics")

# Load project
project_name <- "mouseBrain"

proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)


# 1. Subset the project (already done)
clusters_to_remove <- c("C1", "C2", "C3", "C4", "C20")
cells_to_keep <- proj_ALL$cellNames[!proj_ALL$Clusters_Combined %in% clusters_to_remove]

proj_SUB <- subsetArchRProject(
  ArchRProj = proj_ALL,
  cells = cells_to_keep,
  outputDirectory = paste0(project_name, "_filtered"),
  dropCells = TRUE,
  force = TRUE
)

# Re-run dimensionality reduction and embeddings on filtered data
proj_SUB <- addCombinedDims(proj_SUB, reducedDims = c("LSI_ATAC", "LSI_RNA"), name = "LSI_Combined")
proj_SUB <- addUMAP(proj_SUB, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj_SUB <- addUMAP(proj_SUB, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj_SUB <- addUMAP(proj_SUB, reducedDims = "LSI_Combined", dimsToUse = 1:65, name = "UMAP_Combined", minDist = 0.9, force = TRUE)
# Re-run clustering
proj_SUB <- addClusters(proj_SUB, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 1.0, force = TRUE)
proj_SUB <- addClusters(proj_SUB, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 1.0, force = TRUE)
proj_SUB <- addClusters(proj_SUB, reducedDims = "LSI_Combined", dimsToUse = 1:20, name = "Clusters_Combined", resolution = 3.0, force = TRUE)




##Plot and save 
p <- plotEmbedding(
  ArchRProj = proj_SUB,
  colorBy = "cellColData",
  name = "Clusters_Combined",
  embedding = "UMAP_Combined"
)

figure_name <- paste0(project_name, "_CellTypeUMAP_Filtered.png")

ggsave(
  filename = figure_name,
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)


saveArchRProject(
  ArchRProj = proj_SUB,
  outputDirectory = "mouseBrain_filtered",
  load = FALSE
)
