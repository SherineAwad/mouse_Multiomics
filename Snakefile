     
rule all:
         input:
            expand("mouseBrain/ArrowFiles/{Control}.arrow", Control= config['Control']), 
            expand("mouseBrain/ArrowFiles/{KO}.arrow", KO= config['KO']),        
            "mouseBrain_preFilterQC.pdf",
            "mouseBrain_postFilterQC.pdf",
            "mouseBrain_SamplesUMAP_byCluster.png",
            "mouseBrain_ClustersUMAP.png",
            "mouseBrain_perClustersnUMI.png", 
            "mouseBrain/Plots/UMAP_GeneScore_AllMarkers.pdf",
            "Rplots.pdf",
            "mouseBrain_CellTypeUMAP_annotated.png",
            "mouseBrain_CellTypeUMAP_Filtered.png", 
            "mouseBrain_filtered_CellTypeUMAP_annotated.png"

 
rule preprocess: 
        input:
             "TH1_filtered_feature_bc_matrix.h5",
             "TH2_filtered_feature_bc_matrix.h5",  
        output:
           expand("mouseBrain/ArrowFiles/{Control}.arrow", Control= config['Control']),
           expand("mouseBrain/ArrowFiles/{KO}.arrow", KO= config['KO']), 
        shell: 
            """
            Rscript preprocess.R 
            """ 
     
rule prefilter: 
        input: 
           expand("mouseBrain/ArrowFiles/{Control}.arrow", Control= config['Control']),
           expand("mouseBrain/ArrowFiles/{KO}.arrow", KO= config['KO']),
        output: 
           "mouseBrain_preFilterQC.pdf",
           "mouseBrain_postFilterQC.pdf"
        shell: 
           """
           Rscript preFilter.R
           Rscript filter.R 
           """
      
         
rule UMAP: 
     input: 
           "mouseBrain_preFilterQC.pdf",
           "mouseBrain_postFilterQC.pdf" 
     output: 
           "mouseBrain_SamplesUMAP_byCluster.png",
           "mouseBrain_ClustersUMAP.png",
           "mouseBrain_perClustersnUMI.png",
     shell: 
          """
          Rscript addUMAP.R
          """
rule Inspect: 
     input: 
           expand("mouseBrain/ArrowFiles/{Control}.arrow", Control= config['Control']),
           expand("mouseBrain/ArrowFiles/{KO}.arrow", KO= config['KO']),
     output: 
        "mouseBrain/Plots/QC_on_metrics.pdf",
     shell:  
       """
       Rscript inspect.R 
       """


rule callDGE: 
     input:
         "mouseBrain_SamplesUMAP_byCluster.png",
         "mouseBrain_ClustersUMAP.png",
         "mouseBrain_perClustersnUMI.png", 
     output: 
         "mouseBrain/Plots/UMAP_GeneScore_AllMarkers.pdf",
         "Rplots.pdf" 
     shell:
         """
         Rscript callDGE.R
         Rscript plotMarkers.R
         """


rule annotate: 
     input:
        "mouseBrain/Plots/UMAP_GeneScore_AllMarkers.pdf",
     output: 
        "mouseBrain_CellTypeUMAP_annotated.png"
     shell: 
        """
        Rscript annotate.R mouseBrain cluster_annotations.csv 
        """ 

rule remove_outliers: 
    input: 
       "mouseBrain_CellTypeUMAP_annotated.png"
    output: 
       "mouseBrain_CellTypeUMAP_Filtered.png" 
    shell: 
       """ 
       Rscript remove_outlier.R
       """


rule reannotate: 
    input: 
       "mouseBrain_CellTypeUMAP_Filtered.png"
    output: 
      "mouseBrain_filtered_CellTypeUMAP_annotated.png" 
    shell: 
      """ 
      Rscript annotate.R mouseBrain_filtered cluster_reannotations.csv
      """  
