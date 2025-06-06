# Mouse Brain Multi-Omics Analysis

## Overview

This project analyzes **single-cell RNA-seq (scRNA-seq)** and **single-cell ATAC-seq (scATAC-seq)** data from mouse brain tissue using the **[ArchR](https://www.archrproject.com/)** framework.

## Tools

- **ArchR**: For large-scale single-cell ATAC-seq analysis and integration with RNA-seq.
- **Snakemake**: (optional, for reproducibility and workflow management)
- **R** and **Bioconductor**: For scripting and analysis.
- **Scanpy**: (used as reference or for RNA-only processing)

## Data

- Organism: *Mus musculus* (Mouse)
- Tissue: Brain
- Modalities: scRNA-seq and scATAC-seq


## UMAP and Cluster Plots


### nUMI per Cluster
![nUMI per Cluster](mouseBrain_perClustersnUMI.png)


### Clustering Overview
![Clusters UMAP](mouseBrain_ClustersUMAP.png)

### UMAP by Cluster
![UMAP by Cluster](mouseBrain_SamplesUMAP_byCluster.png)

### UMAP by Sample
![UMAP by Sample](mouseBrain_SamplesUMAP_bySample.png)


### Marker Genes UMAPs 
[MARKER GENES UMAP](Rplots.pdf)

### Annotations 
![Annotations](mouseBrain_CellTypeUMAP_annotated.png)

### Filtered: We removed cluster 1,2, 3, 4, and 20  then reClustered
![Filtered](mouseBrain_CellTypeUMAP_Filtered.png)

### Marker Genes UMAPs after removing outlier Clusters
[MARKER GENES UMAP After removing outliers](Rplotsfiltered.pdf)

### Annotations after removing outlier Clusters
![Annotations After removing outliers](mouseBrain_filtered_CellTypeUMAP_annotated.png) 



### Cell Counts Per Sample

| Sample  | Cell Count |
|---------|------------|
| Control | 3686       |
| KO      | 3757       |




