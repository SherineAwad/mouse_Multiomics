     
rule all:
         input:
            expand("mouseBrain/ArrowFiles/{Control}.arrow", Control= config['Control']), 
            expand("mouseBrain/ArrowFiles/{KO}.arrow", KO= config['KO']),        
            "mouseBrain_preFilterQC.pdf",
            "mouseBrain_postFilterQC.pdf",
            "mouseBrain_SamplesUMAP.pdf",
            "mouseBrain_clustersUMAP.pdf", 
            "mouseBrain_perClustersnUMI.pdf",
            "mouseBrain/Plots/QC_on_metrics.pdf", 

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
          "mouseBrain_SamplesUMAP.pdf", 
          "mouseBrain_clustersUMAP.pdf", 
          "mouseBrain_perClustersnUMI.pdf",   
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
        "mouseBrain_SamplesUMAP.pdf",
        "mouseBrain_clustersUMAP.pdf",
        "mouseBrain_perClustersnUMI.pdf",
     output: 
         "UMAP_GeneScore_AllMarkers.pdf"  
     shell:
         """
         Rscript callDGE.R
         """

