library("tidyverse")
library("DESeq2")
library("ggplot2")
library("hexbin")
library("genefilter")
library("pheatmap")
library("clusterProfiler")
library("gghighlight")
library("RColorBrewer")
library("Cairo")

#Gene Ontology (GO) classifies functions along molecular function (molecular activities), cellular component (where gene products are active), and biological process (pathways and larger processes)
#We have 60609 genes detected, with 3962 that were differentially expressed (padj < 0.05) across samples 
#We selected padj < 0.05 as the differential genes 
#Now we need to evaluate whether those 3962 genes are over-represented in the gene set using hypergeometric distribution 
#and use gene set enrichment analysis (GSEA) to find genes where the difference is large 
#Input into clusterProfiler is a ranked list of genes (i.e., ControlvsCRISPRi_results.csv")

#selected 447 genes with padj < 0.05 
WTvex9_genes_padj <- as.data.frame(read.csv("WTvex9.csv"))
WTvex9_genes_padj_entrez <- as.character(WTvex9_genes_padj$entrez)

#selected 1944 genes with padj < 0.05
WTvex20_genes_padj <- as.data.frame(read.csv("WTvex20.csv"))
WTvex20_genes_padj_entrez <- as.character(WTvex20_genes_padj$entrez)

#selected 2388 genes with padj < 0.05
ex9vex20_genes_padj <- as.data.frame(read.csv("ex9vex20.csv"))
ex9vex20_genes_padj_entrez <- as.character(ex9vex20_genes_padj$entrez)


#filter by l2fc
#selected 1697 genes with log2foldchange outside of (+/-)1.5 and padj < 0.05 
WTvex9_genes_log2fc_1.5 <- WTvex9_genes_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
WTvex9_genes_log2fc_1.5_entrez <- as.character(WTvex9_genes_log2fc_1.5$entrez)
WTvex9_Comparison_gene_list_df<- as.data.frame(WTvex9_genes_log2fc_1.5$entrez)
write.csv(WTvex9_genes_log2fc_1.5, file = "WTvex9_Comparison_Sig_entrez_list.csv")

#selected 1370 genes with log2foldchange outside of (+/-)1.5 and padj < 0.05
WTvex20_genes_log2fc_1.5 <- WTvex20_genes_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
WTvex20_genes_log2fc_1.5_entrez <- as.character(WTvex20_genes_log2fc_1.5$entrez)
WTvex20_genes_log2fc_1.5_entrez_df <- as.data.frame(WTvex20_genes_log2fc_1.5$entrez)
write.csv(WTvex20_genes_log2fc_1.5, file = "WTvex20_Comparison_Sig_entrez_list")


#selected 10684 genes with log2foldchange outside of (+/-)1.5 and padj < 0.05
ex9vex20_genes_log2fc_1.5 <- ex9vex20_genes_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
ex9vex20_genes_log2fc_1.5_entrez <- as.character(ex9vex20_genes_log2fc_1.5$entrez)
ex9vex20_genes_log2fc_1.5_entrez_df <- as.data.frame(ex9vex20_genes_log2fc_1.5$entrez)


######Biological Process######
WTvex9_GOBP <- enrichGO(WTvex9_genes_log2fc_1.5_entrez,
                       ont = "BP",
                       OrgDb="org.Hs.eg.db",
                       pvalueCutoff = 0.05,
                       keyType = "ENTREZID",
                       pAdjustMethod = "BH")

dotplot(WTvex9_GOBP,
        title = "WTvex9 Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

WTvex20_GOBP <- enrichGO(WTvex20_genes_log2fc_1.5_entrez,
                        ont = "BP",
                        OrgDb="org.Hs.eg.db",
                        pvalueCutoff = 0.05,
                        keyType = "ENTREZID",
                        pAdjustMethod = "BH")

dotplot(WTvex20_GOBP,
        title = "WTvex20 Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

ex9vex20_GOBP <- enrichGO(ex9vex20_genes_log2fc_1.5_entrez,
                            ont = "BP",
                            OrgDb="org.Hs.eg.db",
                            pvalueCutoff = 0.05,
                            keyType = "ENTREZID",
                            pAdjustMethod = "BH")

dotplot(ex9vex20_GOBP,
        title = "ex9vex20 Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

#GOMF
WTvex9_GOMF <- enrichGO(WTvex9_genes_log2fc_1.5_entrez,
                       ont = "MF",
                       OrgDb="org.Hs.eg.db",
                       pvalueCutoff = 0.05,
                       keyType = "ENTREZID",
                       pAdjustMethod = "BH")

dotplot(WTvex9_GOMF,
        title = "WTvex9 Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

WTvex20_GOMF <- enrichGO(WTvex20_genes_log2fc_1.5_entrez,
                        ont = "MF",
                        OrgDb="org.Hs.eg.db",
                        pvalueCutoff = 0.05,
                        keyType = "ENTREZID",
                        pAdjustMethod = "BH")

dotplot(WTvex20_GOMF,
        title = "WTvex20Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

ex9vex20_GOMF <- enrichGO(ex9vex20_genes_log2fc_1.5_entrez,
                         ont = "MF",
                         OrgDb="org.Hs.eg.db",
                         pvalueCutoff = 0.05,
                         keyType = "ENTREZID",
                         pAdjustMethod = "BH")

dotplot(ex9vex20_GOMF,
        title = "ex9vex20 Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

###Enrich KEGG
WTvex9_KEGG <- enrichKEGG(WTvex9_genes_log2fc_1.5_entrez,
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = .05,
                        pAdjustMethod = "BH")

dotplot(WTvex9_KEGG,
        title = "WTvex9_KEGG Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 5)

WTvex20_KEGG <- enrichKEGG(WTvex20_genes_log2fc_1.5_entrez,
                          organism = "hsa",
                          keyType = "kegg",
                          pvalueCutoff = .05,
                          pAdjustMethod = "BH")

dotplot(WTvex20_KEGG,
        title = "WTvex20_KEGG Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

ex9vex20_KEGG <- enrichKEGG(ex9vex20_genes_log2fc_1.5_entrez,
                          organism = "hsa",
                          keyType = "kegg",
                          pvalueCutoff = .05,
                          pAdjustMethod = "BH")

dotplot(ex9vex20_KEGG,
        title = "ex9vex20_KEGG Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)