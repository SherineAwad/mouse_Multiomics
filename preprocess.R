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

setwd("/nfs/turbo/umms-thahoang/sherine/mouse_Multiomics")
addArchRThreads(threads = 4) 
addArchRGenome("mm10")
project_name ="mouseBrain"
atacFiles <- c("Control" = "TH1_atac_fragments.tsv.gz", "KO" = "TH2_atac_fragments.tsv.gz")


rnaFiles <- c("Control" = "TH1_filtered_feature_bc_matrix.h5", "KO" = "TH2_filtered_feature_bc_matrix.h5")


names(atacFiles)
names(rnaFiles)
all.equal(names(atacFiles), names(rnaFiles))
ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE)

ArrowFiles <- c("Control.arrow","KO.arrow")
project_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "mouseBrain", copyArrows = FALSE)
#need to check indeces of ArrowFiles
proj_A <- ArchRProject(ArrowFiles[1], outputDirectory = "Control", copyArrows = FALSE)
proj_B <- ArchRProject(ArrowFiles[2], outputDirectory = "KO", copyArrows = FALSE)

seRNA_A <- import10xFeatureMatrix(input="TH1_filtered_feature_bc_matrix.h5", names="Control")
seRNA_B <- import10xFeatureMatrix(input="TH2_filtered_feature_bc_matrix.h5", names="KO")


##Check seqlevels(seRNA_A) and it shows the following lines are not needed 
#gr <- rowRanges(seRNA_A)
#seqlevels(gr) <- paste0("chr", seqlevels(gr))
#rowRanges(seRNA_A) <- gr
k <- which(seqnames(rowRanges(seRNA_A)) %in% seqnames(proj_A@genomeAnnotation$chromSizes) == T)
seRNA_A = seRNA_A[k,]
k = which(rownames(proj_A@cellColData) %in% colnames(seRNA_A) == T)
proj_A = proj_A[k]

##Check seqlevels(seRNA_B) and it shows the following lines are not needed 
#gr <- rowRanges(seRNA_B)
#seqlevels(gr) <- paste0("chr", seqlevels(gr))
#rowRanges(seRNA_B) <- gr
k <- which(seqnames(rowRanges(seRNA_B)) %in% seqnames(proj_B@genomeAnnotation$chromSizes) == T)
seRNA_B = seRNA_B[k,]
k = which(rownames(proj_B@cellColData) %in% colnames(seRNA_B) == T)
proj_B = proj_B[k]



seRNAcombined <- cbind(assay(seRNA_A), assay(seRNA_B))
seRNA_all <- SummarizedExperiment(assays = list(counts = seRNAcombined), rowRanges = rowRanges(seRNA_A))
project_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "mouseBrain", copyArrows = FALSE)
proj_ALL <- addGeneExpressionMatrix(input = project_ALL, seRNA = seRNA_all)
saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "mouseBrain", load = FALSE)






