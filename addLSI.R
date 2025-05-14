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

