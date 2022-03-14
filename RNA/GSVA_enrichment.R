#attempt at using GSVA enrichment on my RNA-seq data
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GSVA")
library(GSVA)
library(GSEABase)
#Using the normalized reads subset from Total_DESeq.R in the clustering portion
#all diff expressed genes from the three comparisons across genotypes
norm_subset_mat <- as.matrix(norm_total_counts_sig_subset)

#GSVA
gsva(norm_subset_mat)
