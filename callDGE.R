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

addArchRThreads(threads = 4)
addArchRGenome("mm10")
setwd("/nfs/turbo/umms-thahoang/sherine/mouse_Multiomics")
project_name ="mouseBrain"
proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

markerGenes <- readLines("markers.txt")


featuresControls <- getMarkerFeatures(
ArchRProj = proj_ALL,
useMatrix = "GeneExpressionMatrix",
groupBy = "Clusters_Combined", bias = c("Gex_nUMI","Gex_nGenes"),
testMethod = "wilcoxon"
)

figure_name = project_name
figure_name <- paste(figure_name,"_ControlsmarkersHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
subsetSE <- featuresControls[which(rowData(featuresControls)$name %in% markerGenes),]
ControlsmarkersHeatmap <- plotMarkerHeatmap(seMarker = subsetSE,cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5",plotLog2FC = TRUE)
draw(ControlsmarkersHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

# Create a list to store all the plots
allPlots <- list()

# Get available gene names
availableGenes <- getFeatures(proj_ALL, useMatrix = "GeneScoreMatrix", select = NULL)

# Loop through genes and generate plots
for (gene in markerGenes) {
  if (gene %in% availableGenes) {
    p <- plotEmbedding(
      ArchRProj = proj_ALL,
      colorBy = "GeneScoreMatrix",
      name = gene,
      embedding = "UMAP_Combined",
    )
    allPlots[[gene]] <- p
  } else {
    message(paste("Skipping gene:", gene, "- not found in GeneScoreMatrix"))
  }
}

# Save all plots in one PDF
plotPDF(allPlots, name = "UMAP_GeneScore_AllMarkers.pdf", ArchRProj = proj_ALL, addDOC = FALSE, width = 5, height = 5)


# Get available gene names from the GeneScoreMatrix
availableGenes <- getFeatures(proj_ALL, useMatrix = "GeneScoreMatrix", select = NULL)

# Loop through genes and plot only if present
for (gene in markerGenes) {
  if (gene %in% availableGenes) {
    p <- plotEmbedding(
      ArchRProj = proj_ALL,
      colorBy = "GeneScoreMatrix",
      name = gene,
      embedding = "UMAP_Combined",
    )
    plotPDF(p, name = paste0("UMAP_GeneScore_", gene, ".pdf"), ArchRProj = proj_ALL, addDOC = FALSE)
  } else {
    message(paste("Skipping gene:", gene, "- not found in GeneScoreMatrix"))
  }
}


saveArchRProject(ArchRProj = proj_ALL, outputDirectory = project_name, load = FALSE)
