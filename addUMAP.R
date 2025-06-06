library(ArchR)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(pheatmap)
library(chromVARmotifs)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(irlba)


addArchRThreads(threads = 4) 
addArchRGenome("mm10")

setwd("/nfs/turbo/umms-thahoang/sherine/mouse_Multiomics")
project_name ="mouseBrain"

proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

#LSI-ATAC
proj_ALL <- addIterativeLSI(
    ArchRProj = proj_ALL,
    clusterParams = list(
      resolution = 0.2,
      sampleCellsPre = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix",
    depthCol = "nFrags",
    name = "LSI_ATAC"
)


#LSI-RNA
proj_ALL <- addIterativeLSI(
    ArchRProj = proj_ALL,
    clusterParams = list(
      resolution = 0.2,
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix",
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA"
)

#-----------------------------------
proj_ALL <- addCombinedDims(proj_ALL, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined", force =TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_Combined", dimsToUse = 1:55, name = "UMAP_Combined", minDist = 0.8, force = TRUE)
#-----------------------------
#Add Clusters
#----------------------------
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 1.0, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 1.0, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_Combined",dimsToUse = 1:20, name = "Clusters_Combined", resolution = 3.0, force = TRUE)

# Violin Plot: perClusters nUMI
figure_name <- paste0(project_name, "_perClustersnUMI.png")
png(filename = figure_name, width = 1200, height = 800, res = 150)
p <- plotGroups(
  ArchRProj = proj_ALL,
  groupBy = "Clusters_Combined",
  colorBy = "cellColData",
  name = "Gex_nUMI",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
print(p)
dev.off()

# UMAP by Sample and by Clusters (each saved separately)
figure_name1 <- paste0(project_name, "_SamplesUMAP_bySample.png")
png(filename = figure_name1, width = 1200, height = 800, res = 150)
p1 <- plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Combined"
)
print(p1)
dev.off()

figure_name2 <- paste0(project_name, "_SamplesUMAP_byCluster.png")
png(filename = figure_name2, width = 1200, height = 800, res = 150)
p2 <- plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "Clusters_Combined",
  embedding = "UMAP_Combined"
)
print(p2)
dev.off()

# UMAP plots for Clusters: ATAC, RNA, Combined in a single PNG
figure_name3 <- paste0(project_name, "_ClustersUMAP.png")
png(filename = figure_name3, width = 1800, height = 600, res = 150)

p1 <- plotEmbedding(proj_ALL, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors = FALSE, labelMeans = FALSE)
p2 <- plotEmbedding(proj_ALL, name = "Clusters_RNA", embedding = "UMAP_RNA", size = 1.5, labelAsFactors = FALSE, labelMeans = FALSE)
p3 <- plotEmbedding(proj_ALL, name = "Clusters_Combined", embedding = "UMAP_Combined", size = 1.5, labelAsFactors = FALSE, labelMeans = FALSE)

# Apply consistent theme to all
p_customized <- lapply(list(p1, p2, p3), function(x) {
  x +
    theme_ArchR(baseSize = 6.5) +
    theme(
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    guides(color = guide_none(), fill = guide_none())
})

# Combine and print
final_plot <- do.call(cowplot::plot_grid, c(list(ncol = 3), p_customized))
print(final_plot)
dev.off()

# Save the project
saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "mouseBrain", load = FALSE)




