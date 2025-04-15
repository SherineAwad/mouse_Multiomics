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




# Filter the cellColData to remove rows with NA values in the specified columns
cleaned_cellColData <- subset(proj_ALL@cellColData, 
                              !is.na(TSSEnrichment) & 
                              !is.na(log10(nFrags)) & 
                              !is.na(Gex_nUMI) & 
                              !is.na(Gex_nGenes))

# Now, create a new ArchR project with the cleaned cellColData
proj_ALL_clean <- proj_ALL
proj_ALL_clean@cellColData <- cleaned_cellColData

# Proceed with your plotting
figure_name <- paste(project_name, "_preFilterQC.pdf", sep="")
pdf(file = figure_name, width = 12, height = 8)

# Plot 1: TSSEnrichment (violin + boxplot)
p1 <- plotGroups(
    ArchRProj = proj_ALL_clean,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)

# Plot 2: log10(nFrags) (violin + boxplot)
p2 <- plotGroups(
    ArchRProj = proj_ALL_clean,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)

# Plot 3: Gex_nUMI (violin + boxplot)
p3 <- plotGroups(
    ArchRProj = proj_ALL_clean,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "Gex_nUMI",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)

# Plot 4: Gex_nGenes (violin + boxplot)
p4 <- plotGroups(
    ArchRProj = proj_ALL_clean,
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

saveArchRProject(ArchRProj = proj_ALL_clean, outputDirectory = "mouseBrain", load = FALSE)


