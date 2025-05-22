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

# Read marker genes
markerGenes <- readLines("markers.txt")
markerGenes <- trimws(markerGenes)

plot_gene_umap_zscore <- function(proj, gene_symbol) {
  message("Plotting gene: ", gene_symbol)
  flush.console()
  
  # Get expression matrix
  expr_matrix <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
  log_norm_expr <- assay(expr_matrix)
  gene_names <- rowData(expr_matrix)$name
  rownames(log_norm_expr) <- gene_names
  
  # Find gene index case-insensitive
  gene_idx <- which(tolower(rownames(log_norm_expr)) == tolower(gene_symbol))
  if(length(gene_idx) == 0) {
    warning(paste("Gene", gene_symbol, "not found in expression matrix. Skipping."))
    return(NULL)
  }
  
  expr_values <- log_norm_expr[gene_idx, ]
  
  # Get embedding
  embedding <- getEmbedding(proj, embedding = "UMAP_Combined", returnDF = TRUE)
  
  # Check embedding column names and rename for plotting
  print(colnames(embedding))
  if(ncol(embedding) == 2) {
    colnames(embedding) <- c("UMAP_1", "UMAP_2")
  } else {
    stop("Embedding does not have exactly 2 columns! Please check embedding name or data.")
  }
  
  # Make sure cell names match between embedding and expression matrix
  common_cells <- intersect(rownames(embedding), colnames(log_norm_expr))
  if(length(common_cells) == 0) stop("No common cells between embedding and expression matrix!")
  
  # Subset to common cells
  embedding_sub <- embedding[common_cells, , drop = FALSE]
  expr_sub <- expr_values[common_cells]
  
  # Z-score normalization of expression
  expr_zscore <- (expr_sub - mean(expr_sub)) / sd(expr_sub)
  
  # Data frame for plotting
  df <- data.frame(embedding_sub, ExpressionZ = expr_zscore)
  # Plot UMAP colored by z-scored expression
p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = ExpressionZ)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(
    colors = c("#add8e6", "#8B0000"),  # Deep red to light blue
    limits = c(min(df$ExpressionZ, na.rm = TRUE), max(df$ExpressionZ, na.rm = TRUE))
  ) +
  theme_minimal() +
  ggtitle(paste0("UMAP: Z-score normalized ", gene_symbol))

print(p)

# Save PDF using ArchR helper
plotPDF(
  p,
  name = paste0("UMAP_Zscore_", gene_symbol, ".pdf"),
  ArchRProj = proj,
  addDOC = FALSE,
  width = 5,
  height = 5
)

 
 



}

# Loop through marker genes and plot
for (gene in markerGenes) {
  plot_gene_umap_zscore(proj_ALL, gene)
}

