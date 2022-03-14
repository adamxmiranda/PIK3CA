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

###Comparisons for doing shared pathway analyses
total_treatment_v2 <- c("ex9_MCF10A", "ex9_MCF10A", "ex9_MCF10A",
                        "WT_MCF10A", "WT_MCF10A", "WT_MCF10A",
                        "WT_MCF7", "WT_MCF7", "WT_MCF7",
                        "ex9_MCF7", "ex9_MCF7", "ex9_MCF7",
                        "ex20_MCF10A", "ex20_MCF10A", "ex20_MCF10A",
                        "WT_HTert", "WT_HTert", "WT_HTert",
                        "ex9_HTert", "ex9_HTert", "ex9_HTert",
                        "ex20_HTert", "ex20_HTert", "ex20_HTert",
                        "ex20_T47D", "ex20_T47D", "ex20_T47D")
total_coldata_v2 <- data.frame(total_treatment_v2)
total_dds_cellSpec <- DESeqDataSetFromMatrix(countData= total_counts,
                                             colData= total_coldata_v2,
                                             design= ~ total_treatment_v2) 
total_dds_cellSpec <- DESeq(total_dds_cellSpec)
resultsNames(total_dds_cellSpec)

#Mutant v Mutants
HTert_ex9vex20_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "ex9_HTert", "ex20_HTert"))
cancer_ex9vex20_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "ex9_MCF7", "ex20_T47D"))
MCF10A_ex9vex20_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "ex9_MCF10A", "ex20_MCF10A"))

#Mutant v WT
HTert_WTvex9_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "WT_HTert", "ex9_HTert"))
cancer_WTvex9_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "WT_MCF7", "ex9_MCF7"))
MCF10A_WTvex9_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "WT_MCF10A", "ex9_MCF10A"))
HTert_WTvex20_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "WT_HTert", "ex20_HTert"))
cancer_WTvex20_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "WT_MCF7", "ex20_T47D"))
MCF10A_WTvex20_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "WT_MCF10A", "ex20_MCF10A"))

#Mutant v Mutants
HTert_ex9vex20_res_df <- as.data.frame(HTert_ex9vex20_res)
HTert_ex9vex20_res_df_filt <- HTert_ex9vex20_res_df %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>%
  filter(padj < 0.05)
cancer_ex9vex20_res_df <- as.data.frame(cancer_ex9vex20_res)
cancer_ex9vex20_res_df_filt <- cancer_ex9vex20_res_df %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>%
  filter(padj < 0.05)
MCF10A_ex9vex20_res_df <- as.data.frame(MCF10A_ex9vex20_res)
MCF10A_ex9vex20_res_df_filt <- MCF10A_ex9vex20_res_df %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>%
  filter(padj < 0.05)

#Mutant v WT
HTert_WTvex9_res_df <- as.data.frame(HTert_WTvex9_res)
HTert_WTvex9_res_df_filt <- HTert_WTvex9_res_df %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>%
  filter(padj < 0.05)
cancer_WTvex9_res_df <- as.data.frame(cancer_WTvex9_res)
cancer_WTvex9_res_df_filt <- cancer_WTvex9_res_df %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>%
  filter(padj < 0.05)
MCF10A_WTvex9_res_df <- as.data.frame(MCF10A_WTvex9_res)
MCF10A_WTvex9_res_df_filt <- MCF10A_WTvex9_res_df %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>%
  filter(padj < 0.05)

HTert_WTvex20_res_df <- as.data.frame(HTert_WTvex20_res)
HTert_WTvex20_res_df_filt <- HTert_WTvex20_res_df %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>%
  filter(padj < 0.05)
cancer_WTvex20_res_df <- as.data.frame(cancer_WTvex20_res)
cancer_WTvex20_res_df_filt <- cancer_WTvex20_res_df %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>%
  filter(padj < 0.05)
MCF10A_WTvex20_res_df <- as.data.frame(MCF10A_WTvex20_res)
MCF10A_WTvex20_res_df_filt <- MCF10A_WTvex20_res_df %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5) %>%
  filter(padj < 0.05)

#Get Entrez IDs
library(org.Hs.eg.db)
gene_lists <- list(HTert_ex9vex20_res_df_filt,
                   cancer_ex9vex20_res_df_filt,
                   MCF10A_ex9vex20_res_df_filt)
for (i in 1:length(gene_lists)) {
  gene_list <- rownames(gene_lists[[i]])
  ens.str <- substr(x = gene_list, 1 , 15)
  entrez_col <- mapIds(org.Hs.eg.db,
                       keys = ens.str,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals= "first")
  gene_lists[[i]] <- cbind(gene_lists[[i]], entrez_col)
}
#WTvMut
library(org.Hs.eg.db)
WTvMut_gene_lists <- list(HTert_WTvex9_res_df_filt,
                          cancer_WTvex9_res_df_filt,
                          MCF10A_WTvex9_res_df_filt,
                          HTert_WTvex20_res_df_filt,
                          cancer_WTvex20_res_df_filt,
                          MCF10A_WTvex20_res_df_filt)
for (i in 1:length(WTvMut_gene_lists)) {
  gene_list <- rownames(WTvMut_gene_lists[[i]])
  ens.str <- substr(x = gene_list, 1 , 15)
  entrez_col <- mapIds(org.Hs.eg.db,
                       keys = ens.str,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals= "first")
  WTvMut_gene_lists[[i]] <- cbind(WTvMut_gene_lists[[i]], entrez_col)
}
WTvMut_goodNames <- c("HTert_WTvex9_res_df_filt",
                      "cancer_WTvex9_res_df_filt",
                      "MCF10A_WTvex9_res_df_filt",
                      "HTert_WTvex20_res_df_filt",
                      "cancer_WTvex20_res_df_filt",
                      "MCF10A_WTvex20_res_df_filt")
names(WTvMut_gene_lists) <- WTvMut_goodNames


goodNames <- c("HTert_ex9vex20_res_df_filt",
               "cancer_ex9vex20_res_df_filt",
               "MCF10A_ex9vex20_res_df_filt")
names(gene_lists) <- goodNames
library(clusterProfiler)
HTert_ex9vex20_entrez <- gene_lists$HTert_ex9vex20_res_df_filt$entrez_col
HTert_ex9vex20_entrez <- na.omit(HTert_ex9vex20_entrez)
HTert_ex9vex20_KEGG <- enrichKEGG(gene = HTert_ex9vex20_entrez,
                                  organism = "hsa",
                                  keyType = "kegg",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05)
cancer_ex9vex20_entrez <- gene_lists$cancer_ex9vex20_res_df_filt$entrez_col
cancer_ex9vex20_entrez <- na.omit(cancer_ex9vex20_entrez)
cancer_ex9vex20_KEGG <- enrichKEGG(gene = cancer_ex9vex20_entrez,
                                   organism = "hsa",
                                   keyType = "kegg",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05)
MCF10A_ex9vex20_entrez <- gene_lists$MCF10A_ex9vex20_res_df_filt$entrez_col
MCF10A_ex9vex20_entrez <- na.omit(MCF10A_ex9vex20_entrez)
MCF10A_ex9vex20_KEGG <- enrichKEGG(gene = MCF10A_ex9vex20_entrez,
                                   organism = "hsa",
                                   keyType = "kegg",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05)

#WTvMut
HTert_WTvex9_entrez <- WTvMut_gene_lists$HTert_WTvex9_res_df_filt$entrez_col
HTert_WTvex9_entrez <- na.omit(HTert_WTvex9_entrez)
HTert_WTvex9_KEGG <- enrichKEGG(gene = HTert_WTvex9_entrez,
                                  organism = "hsa",
                                  keyType = "kegg",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05)
cancer_WTvex9_entrez <- WTvMut_gene_lists$cancer_WTvex9_res_df_filt$entrez_col
cancer_WTvex9_entrez <- na.omit(cancer_WTvex9_entrez)
cancer_WTvex9_KEGG <- enrichKEGG(gene = cancer_WTvex9_entrez,
                                   organism = "hsa",
                                   keyType = "kegg",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05)
MCF10A_WTvex9_entrez <- WTvMut_gene_lists$MCF10A_WTvex9_res_df_filt$entrez_col
MCF10A_WTvex9_entrez <- na.omit(MCF10A_WTvex9_entrez)
MCF10A_WTvex9_KEGG <- enrichKEGG(gene = MCF10A_WTvex9_entrez,
                                   organism = "hsa",
                                   keyType = "kegg",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05)
HTert_WTvex20_entrez <- WTvMut_gene_lists$HTert_WTvex20_res_df_filt$entrez_col
HTert_WTvex20_entrez <- na.omit(HTert_WTvex20_entrez)
HTert_WTvex20_KEGG <- enrichKEGG(gene = HTert_WTvex20_entrez,
                                  organism = "hsa",
                                  keyType = "kegg",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05)
cancer_WTvex20_entrez <- WTvMut_gene_lists$cancer_WTvex20_res_df_filt$entrez_col
cancer_WTvex20_entrez <- na.omit(cancer_WTvex20_entrez)
cancer_WTvex20_KEGG <- enrichKEGG(gene = cancer_WTvex20_entrez,
                                   organism = "hsa",
                                   keyType = "kegg",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05)
MCF10A_WTvex20_entrez <- WTvMut_gene_lists$MCF10A_WTvex20_res_df_filt$entrez_col
MCF10A_WTvex20_entrez <- na.omit(MCF10A_WTvex20_entrez)
MCF10A_WTvex20_KEGG <- enrichKEGG(gene = MCF10A_WTvex20_entrez,
                                   organism = "hsa",
                                   keyType = "kegg",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05)

dotplot(HTert_ex9vex20_KEGG,
        showCategory = 20,
        title = "HTert L2fc(+/-1.5) KEGG")
dotplot(cancer_ex9vex20_KEGG,
        showCategory = 20,
        title = "cancer L2fc(+/-1.5) KEGG")
dotplot(MCF10A_ex9vex20_KEGG,
        showCategory = 20,
        title = "MCF10A L2fc(+/-1.5) KEGG")

dotplot(HTert_WTvex9_KEGG,
        showCategory = 20,
        title = "HTert WTvex9 L2fc(+/-1.5) KEGG")
dotplot(cancer_WTvex9_KEGG,
        showCategory = 20,
        title = "cancer WTvex9 L2fc(+/-1.5) KEGG")
dotplot(MCF10A_WTvex9_KEGG,
        showCategory = 20,
        title = "MCF10A WTvex9 L2fc(+/-1.5) KEGG")

dotplot(HTert_WTvex20_KEGG,
        showCategory = 20,
        title = "HTert WTvex20 L2fc(+/-1.5) KEGG")
dotplot(cancer_WTvex20_KEGG,
        showCategory = 20,
        title = "cancer WTvex20 L2fc(+/-1.5) KEGG")
dotplot(MCF10A_WTvex20_KEGG,
        showCategory = 20,
        title = "MCF10A WTvex20 L2fc(+/-1.5) KEGG")

WTvex9_KEGG_lists <- list(HTert_WTvex9_KEGG,
                          cancer_WTvex9_KEGG,
                          MCF10A_WTvex9_KEGG)
names(WTvex9_KEGG_lists) <- c("HTert_WTvex9_KEGG",
                              "cancer_WTvex9_KEGG",
                              "MCF10A_WTvex9_KEGG")
WTvex20_KEGG_lists <- list(HTert_WTvex20_KEGG,
                          cancer_WTvex20_KEGG,
                          MCF10A_WTvex20_KEGG)
names(WTvex20_KEGG_lists) <- c("HTert_WTvex20_KEGG",
                           "cancer_WTvex20_KEGG",
                           "MCF10A_WTvex20_KEGG")
WTvex9_KEGG_list_df <-lapply(WTvex9_KEGG_lists, as.data.frame)
WTvex20_KEGG_list_df <- lapply(WTvex20_KEGG_lists, as.data.frame)

HTert_ex9vex20_KEGG_df <- as.data.frame(HTert_ex9vex20_KEGG)
cancer_ex9vex20_KEGG_df <- as.data.frame(cancer_ex9vex20_KEGG)
MCF10A_ex9vex20_KEGG_df <- as.data.frame(MCF10A_ex9vex20_KEGG)

sharedPaths <- intersect(intersect(HTert_ex9vex20_KEGG_df$ID,
                                   cancer_ex9vex20_KEGG_df$ID),
                         MCF10A_ex9vex20_KEGG_df$ID)

WTvex9_sharedPaths <- intersect(intersect(WTvex9_KEGG_list_df$HTert_WTvex9_KEGG$ID,
                                          WTvex9_KEGG_list_df$cancer_WTvex9_KEGG$ID),
                                WTvex9_KEGG_list_df$MCF10A_WTvex9_KEGG$ID)
WTvex20_sharedPaths <- intersect(intersect(WTvex20_KEGG_list_df$HTert_WTvex20_KEGG$ID,
                                          WTvex20_KEGG_list_df$cancer_WTvex20_KEGG$ID),
                                WTvex20_KEGG_list_df$MCF10A_WTvex20_KEGG$ID)
WTvex9_sharedPaths
WTvex20_sharedPaths

sharedPaths
sharedPaths_df <- HTert_ex9vex20_KEGG_df %>%
  subset(ID %in% sharedPaths) %>%
  select(ID, Description)
sharedPaths_df

WTvex20_sharedPaths
WTvex20_sharedPaths_df <- WTvex20_KEGG_list_df$HTert_WTvex20_KEGG %>%
  subset(ID %in% WTvex20_sharedPaths) %>%
  select(ID, Description)
WTvex20_sharedPaths