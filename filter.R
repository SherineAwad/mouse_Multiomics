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

proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)


head(proj_ALL@cellColData)

# Match ATAC quality thresholds — up to you if you keep as-is
proj_ALL <- proj_ALL[proj_ALL$TSSEnrichment > 10 &
                     proj_ALL$nFrags > 5000 &
                     !is.na(proj_ALL$Gex_nUMI)]

# Match Scanpy-style RNA QC
proj_ALL <- proj_ALL[proj_ALL$Gex_nGenes > 1000 &
                     proj_ALL$Gex_nGenes < 7000 &
                     proj_ALL$Gex_nUMI > 1500 &
                     proj_ALL$Gex_nUMI < 30000]




#After Filtering 
table(proj_ALL$Sample)


figure_name <- paste(project_name, "_postFilterQC.pdf", sep="")
pdf(file = figure_name, width = 12, height = 8)

# Plot 1: TSSEnrichment (violin + boxplot)
p1 <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)

# Plot 2: log10(nFrags) (violin + boxplot)
p2 <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)

# Plot 3: Gex_nUMI (violin + boxplot)
p3 <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "Gex_nUMI",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)

# Plot 4: Gex_nGenes (violin + boxplot)
p4 <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "Gex_nGenes",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)

# Combine all the plots into a single figure using patchwork
(p1 | p2) / (p3 | p4)

dev.off()




saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "mouseBrain", load = FALSE)







