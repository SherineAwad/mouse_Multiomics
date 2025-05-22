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

# ================================
# Annotate Clusters with Cell Types
# ================================

# Define your cluster-to-celltype mapping here:
# Example (replace with your actual clusters and labels)
cluster_to_celltype <- c(
  "C5" = "Astro",
  "C6" = "Astro",
  "C7" = "Astro",
  "C24" = "Pro_NSC",
  "C25" = "Pro_NSC",
  "C17" = "Oligos", 
  "C18" = "Oligos", 
  "C19" = "Oligos",  
  "C21" = "Neurons",
  "C22" = "Neurons",
  "C23" = "Neurons", 
  "C9" = "KO_Astro",  
  "C10" = "KO_Astro", 
  "C11" = "KO_Astro", 
  "C14" = "NSC", 
  "C12" = "Dediff_Astro",
  "C13" = "Dediff_Astro",
  "C15" = "Dediff_Astro",
  "C16" = "Dediff_Astro" 
  ) 

# Add to ArchR project metadata
proj_ALL$CellType <- cluster_to_celltype[proj_ALL$Clusters_Combined]


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

