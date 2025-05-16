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

p <-  plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "Clusters_Combined",
  embedding = "UMAP_Combined"
)
figure_name = "UMAP_Combined_Clusters.pdf"
plotPDF(p, name = figure_name, ArchRProj = proj_ALL, addDOC = FALSE)

p1 <- plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "TSSEnrichment",
  embedding = "UMAP_Combined"
)

p2 <- plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "nFrags",
  embedding = "UMAP_Combined"
)

plotPDF(p1, p2, name = "QC_on_UMAP.pdf", ArchRProj = proj_ALL, addDOC = FALSE)




p1 <- plotGroups(
  ArchRProj = proj_ALL,
  groupBy = "Clusters_Combined",  
  colorBy = "cellColData",
  name = "TSSEnrichment"
)


p2 <-  plotGroups(
  ArchRProj = proj_ALL,
  groupBy = "Clusters_Combined",
  colorBy = "cellColData",
  name = "nFrags"
)


p3 <- plotGroups(
  ArchRProj = proj_ALL,
  groupBy = "Clusters_Combined",
  colorBy = "cellColData",
  name = "Gex_nUMI"
)

p4 <- plotGroups(
  ArchRProj = proj_ALL,
  groupBy = "Clusters_Combined",
  colorBy = "cellColData",
  name = "Gex_nGenes"
)



p1 <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 <- p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plotPDF(p1, p2,p3, p4, name = "QC_on_metrics.pdf", ArchRProj = proj_ALL, addDOC = FALSE)
dev.off()

