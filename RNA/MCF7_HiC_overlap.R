#Hi-C to gene expression overlap
#only MCF7 data produced thus far
library(tidyverse)
library("DESeq2")
library("ggplot2")
library("hexbin")
library("apeglm")
library("genefilter")
library("pheatmap")
library("cluster")
library("factoextra")
library("AnnotationDbi")
library("org.Hs.eg.db")
library(ggvenn)
#read in the genes enriched in the Hi-C data
MCF7_parental_hic_enriched <- read_tsv(file = "~/Desktop/Total_RNA/parental_enrich_genes.tsv")
MCF7_corrected_hic_enriched <- read_tsv(file = "~/Desktop/Total_RNA/corrected_enrich_genes.tsv")

#objects containing diff expressed genes from this comparison
#from Under_Over.R
cancer_WTvex9_under
cancer_WTvex9_over

#get gene symbols from the diff expressed genes
#first remove the transcript IDs
cancer_WTvex9_under_clean <- tools::file_path_sans_ext(
  c(cancer_WTvex9_under$rowname))
cancer_WTvex9_under_symbols <- mapIds(org.Hs.eg.db,
                                      keys=cancer_WTvex9_under_clean,
                                      column="SYMBOL",
                                      keytype="ENSEMBL",
                                      multiVals="first")
cancer_WTvex9_under_symbols <- cancer_WTvex9_under_symbols[!is.na(cancer_WTvex9_under_symbols)]

cancer_WTvex9_over_clean <- tools::file_path_sans_ext(
  c(cancer_WTvex9_over$rowname))
cancer_WTvex9_over_symbols <- mapIds(org.Hs.eg.db,
                                      keys=cancer_WTvex9_over_clean,
                                      column="SYMBOL",
                                      keytype="ENSEMBL",
                                      multiVals="first")
cancer_WTvex9_over_symbols <- cancer_WTvex9_over_symbols[!is.na(cancer_WTvex9_over_symbols)]

#make venn diagrams on these lists
HiC_overlap_venns <-list("MCF7 E545K expressed" = c(cancer_WTvex9_under_symbols),
                         "MCF7 WT expressed" = c(cancer_WTvex9_over_symbols),
                         "MCF7 E545K interacting" = c(MCF7_parental_hic_enriched$Ensembl.Genes.Desc),
                         "MCF7 WT interacting" = c(MCF7_corrected_hic_enriched$Ensembl.Genes.Desc))
p <- ggvenn(HiC_overlap_venns, c("MCF7 E545K expressed",
                         "MCF7 WT expressed",
                         "MCF7 E545K interacting",
                         "MCF7 WT interacting"),
            stroke_size = 0.5, set_name_size = 2)
p
#get the overlapping gene lists
E545K_expressed_E545K_interacting_overlap <- intersect(cancer_WTvex9_under_symbols, MCF7_parental_hic_enriched$Ensembl.Genes.Desc)
WT_expressed_WT_interacting_overlap <- intersect(cancer_WTvex9_over_symbols, MCF7_corrected_hic_enriched$Ensembl.Genes.Desc)
WT_expressed_E545K_interacting_overlap <- intersect(cancer_WTvex9_over_symbols, MCF7_parental_hic_enriched$Ensembl.Genes.Desc)
E545K_expressed_WT_interacting_overlap <- intersect(cancer_WTvex9_under_symbols, MCF7_corrected_hic_enriched$Ensembl.Genes.Desc)
