# All analysis for the publication: PIK-ing apart hotspots: A functional genomics process for systematic dissection and mutation-specific target discovery in breast cancer PIK3CA hotspot mutations

*Include full citation when prepared

## Organization Guide

### RNA-seq

All RNA-seq information in the publication folder follows analyses performed in the RNA folder in the main repository

#### Total_RNA_Pathway_Specific.ipynb

contains analyses regarding the identification of Differentially expressed genes and the generation of associated figures.

#### TCGA_BRCA_AREG.ipynb

this notebook contains all of the analyses using the PIK3CA mutant RNA-seq samples from the TCGA and METABRIC database.

#### METABRIC.ipynb

This notebood contains additional analysis of the METABRIC RNA data.

### ATAC-seq

All ATAC-seq information in the publication folder follows analyses performed in the ATAC folder in the main repository

#### ATAC-seq.ipynb

This notebook contains all of the main data preparation for the ATAC-seq data.

#### ATAC_Deeptools Heatmaps.ipynb

This notebook contains the analyses that generated the heatmap figure and the analysis of the mutation-preferred clusters from the heatmap.

#### More_ATAC_Figures.ipynb

This notebook contains the code for the generation of MCF-10A cell line specific figures

#### Plotgardener_figure.ipynb

This notebook contains the code for generating the genomic track images of the AREG enhancer region

### CRISPR KO Screen

Gene target selection is found in the DATA_Unite.ipynb notebook.

The majority of these analyses for this experiment are included in the MCF10A_CRISPRscreen.ipynb notebook.

For the evaluation of the top hits in the TCGA data, see the TCGA_BRCA_AREG notebook.

## Other Experimental Data

### Synthetic_Lethal_n_dCas_Assays.ipynb

This notebook contains the analyses of qPCR and the cell viability following perturbation of the AREG enhancer.

*dCas assays were not completed

### Tapestation_suppFig.ipynb

This notebook contains the generation of the enhancer deletion validation supplemental figure.

### antiAREG.ipynb

This notebook contains the analysis of the results from the AREG inhibition experiments.
