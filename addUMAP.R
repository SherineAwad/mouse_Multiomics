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

class(proj_ALL)

names(proj_ALL)

#proj_ALL <- addTileMatrix(proj_ALL, force = TRUE)  # This will recalculate the TileMatrix



#saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "mouseBrain", load = FALSE)



getAvailableMatrices(proj_ALL)



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

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "mouseBrain", load = FALSE)

if (False)
{
#-----------------------------------
proj_ALL <- addCombinedDims(proj_ALL, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_Combined", dimsToUse = 1:20, name = "UMAP_Combined", minDist = 0.4, force = TRUE)

#-----------------------------
#Add Clusters
#----------------------------
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.6, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.6, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_Combined",dimsToUse = 1:20, name = "Clusters_Combined", resolution = 1.2, force = TRUE)

figure_name <- project_name
figure_name <- paste(figure_name,"_perClustersnUMI.pdf", sep="")
pdf(file =figure_name, width=12)
p <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Clusters_Combined",
    colorBy = "cellColData",
    name = "Gex_nUMI",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p
dev.off()

#----------------------------------
figure_name <- project_name
figure_name <- paste(figure_name,"_clustersUMAP.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotEmbedding(proj_ALL, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj_ALL, name = "Clusters_RNA", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj_ALL, name = "Clusters_Combined", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
p1 + p2 + p3 + patchwork::plot_layout(nrow = 1, guides = "collect")
p <- lapply(list(p1,p2,p3), function(x){
  x + guides(color = "none", fill = "none") +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p))
dev.off()

figure_name <- project_name
figure_name <- paste(figure_name,"_SamplesUMAP.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Combined")
p2 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Clusters_Combined", embedding = "UMAP_Combined")
ggAlignPlots(p1, p2, type = "h")
dev.off()

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "mouseBrain", load = FALSE)


}
