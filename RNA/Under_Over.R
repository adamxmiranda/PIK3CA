#looking for over/underexpression patterns in the resSig data
#from Lineage_Clustering.R
#Load the libraries
library("tidyverse")
library("DESeq2")
library("ggplot2")
library("hexbin")
library("apeglm")
library("genefilter")
library("pheatmap")
library("cluster")
library("factoextra")
library("clusterProfiler")
library("gghighlight")
library("RColorBrewer")
library("Cairo")
library(ggvenn)
library(wesanderson)


#subset and overlap over and under expressed genes across genotypic
#comparisons
#WTvex9
#subset overexpressed
cancer_WTvex9_over <- cancer_WTvex9_resSig_padj %>%
  subset(log2FoldChange > 0) %>%
  as.data.frame()
MCF10A_WTvex9_over <- MCF10A_WTvex9_resSig_padj %>%
  subset(log2FoldChange > 0) %>%
  as.data.frame()
HTert_WTvex9_over <- HTert_WTvex9_resSig_padj %>%
  subset(log2FoldChange > 0) %>%
  as.data.frame()
#subset underexpressed
cancer_WTvex9_under <- cancer_WTvex9_resSig_padj %>%
  subset(log2FoldChange < 0) %>%
  as.data.frame()
MCF10A_WTvex9_under <- MCF10A_WTvex9_resSig_padj %>%
  subset(log2FoldChange < 0) %>%
  as.data.frame()
HTert_WTvex9_under <- HTert_WTvex9_resSig_padj %>%
  subset(log2FoldChange < 0) %>%
  as.data.frame()
#create over/under superlist
WTvex9_over_common <- intersect(intersect(cancer_WTvex9_over$rowname,
                                    MCF10A_WTvex9_over$rowname),
                                    HTert_WTvex9_over$rowname)
WTvex9_under_common <- intersect(intersect(cancer_WTvex9_under$rowname,
                                          MCF10A_WTvex9_under$rowname),
                                HTert_WTvex9_under$rowname)
###both are only 2 genes!!!
#WTvex20
#subset overexpressed
cancer_WTvex20_over <- cancer_WTvex20_resSig_padj %>%
  subset(log2FoldChange > 0) %>%
  as.data.frame()
MCF10A_WTvex20_over <- MCF10A_WTvex20_resSig_padj %>%
  subset(log2FoldChange > 0) %>%
  as.data.frame()
HTert_WTvex20_over <- HTert_WTvex20_resSig_padj %>%
  subset(log2FoldChange > 0) %>%
  as.data.frame()
#subset underexpressed
cancer_WTvex20_under <- cancer_WTvex20_resSig_padj %>%
  subset(log2FoldChange < 0) %>%
  as.data.frame()
MCF10A_WTvex20_under <- MCF10A_WTvex20_resSig_padj %>%
  subset(log2FoldChange < 0) %>%
  as.data.frame()
HTert_WTvex20_under <- HTert_WTvex20_resSig_padj %>%
  subset(log2FoldChange < 0) %>%
  as.data.frame()
#create over/under superlist
WTvex20_over_common <- intersect(intersect(cancer_WTvex20_over$rowname,
                                          MCF10A_WTvex20_over$rowname),
                                HTert_WTvex20_over$rowname)
###19 genes
WTvex20_under_common <- intersect(intersect(cancer_WTvex20_under$rowname,
                                           MCF10A_WTvex20_under$rowname),
                                 HTert_WTvex20_under$rowname)
###14 genes
#ex9vex20
#subset overexpressed
cancer_ex9vex20_over <- cancer_ex9vex20_resSig_padj %>%
  subset(log2FoldChange > 0) %>%
  as.data.frame()
MCF10A_ex9vex20_over <- MCF10A_ex9vex20_resSig_padj %>%
  subset(log2FoldChange > 0) %>%
  as.data.frame()
HTert_ex9vex20_over <- HTert_ex9vex20_resSig_padj %>%
  subset(log2FoldChange > 0) %>%
  as.data.frame()
#subset underexpressed
cancer_ex9vex20_under <- cancer_ex9vex20_resSig_padj %>%
  subset(log2FoldChange < 0) %>%
  as.data.frame()
MCF10A_ex9vex20_under <- MCF10A_ex9vex20_resSig_padj %>%
  subset(log2FoldChange < 0) %>%
  as.data.frame()
HTert_ex9vex20_under <- HTert_ex9vex20_resSig_padj %>%
  subset(log2FoldChange < 0) %>%
  as.data.frame()
#create over/under superlist
ex9vex20_over_common <- intersect(intersect(cancer_ex9vex20_over$rowname,
                                           MCF10A_ex9vex20_over$rowname),
                                 HTert_ex9vex20_over$rowname)
###16 genes
ex9vex20_under_common <- intersect(intersect(cancer_ex9vex20_under$rowname,
                                            MCF10A_ex9vex20_under$rowname),
                                  HTert_ex9vex20_under$rowname)
###17 genes

#Let's make some Venns
#using ggvenn
venn_data <- list("WTvex9_over" = c(WTvex9_over_common),
                  "WTvex9_under" = c (WTvex9_under_common),
                  "WTvex20_over" = c(WTvex20_over_common),
                  "WTvex20_under" = c(WTvex20_under_common),
                  "ex9vex20_over" = c(ex9vex20_over_common),
                  "ex9vex20_under" = c(ex9vex20_under_common))
p <- ggvenn(venn_data, c("WTvex9_under",
                         "WTvex20_under",
                         "ex9vex20_over",
                         "ex9vex20_under"),
            stroke_size = 0.5, set_name_size = 2)
p

#Let's try some clusterprofiling
#clean the data of transcript IDs
WTvex9_over_clean <- tools::file_path_sans_ext(
  c(WTvex9_over_common))
WTvex9_under_clean <- tools::file_path_sans_ext(
  c(WTvex9_under_common))
WTvex20_over_clean <- tools::file_path_sans_ext(
  c(WTvex20_over_common))
WTvex20_under_clean <- tools::file_path_sans_ext(
  c(WTvex20_under_common))
ex9vex20_over_clean <- tools::file_path_sans_ext(
  c(ex9vex20_over_common))
ex9vex20_under_clean <- tools::file_path_sans_ext(
  c(ex9vex20_under_common))
#BP
WTvex9_over_common_GOBP <- enrichGO(WTvex9_over_clean,
                                      ont = "BP",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENSEMBL",
                                      pAdjustMethod = "BH")
WTvex9_under_common_GOBP <- enrichGO(WTvex9_under_clean,
                                    ont = "BP",
                                    OrgDb="org.Hs.eg.db",
                                    pvalueCutoff = 0.05,
                                    keyType = "ENSEMBL",
                                    pAdjustMethod = "BH")
WTvex20_over_common_GOBP <- enrichGO(WTvex20_over_clean,
                                     ont = "BP",
                                     OrgDb="org.Hs.eg.db",
                                     pvalueCutoff = 0.05,
                                     keyType = "ENSEMBL",
                                     pAdjustMethod = "BH")
WTvex20_under_common_GOBP <- enrichGO(WTvex20_under_clean,
                                     ont = "BP",
                                     OrgDb="org.Hs.eg.db",
                                     pvalueCutoff = 0.05,
                                     keyType = "ENSEMBL",
                                     pAdjustMethod = "BH")
ex9vex20_over_common_GOBP <- enrichGO(ex9vex20_over_clean,
                                      ont = "BP",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENSEMBL",
                                      pAdjustMethod = "BH")
ex9vex20_under_common_GOBP <- enrichGO(ex9vex20_under_clean,
                                      ont = "BP",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENSEMBL",
                                      pAdjustMethod = "BH")
#BPplots
dotplot(WTvex9_over_common_GOBP,
        title = "WTvex9_over",
        showCategory = 15)
dotplot(WTvex9_under_common_GOBP,
        title = "WTvex9_under",
        showCategory = 15)
dotplot(WTvex20_over_common_GOBP,
        title = "WTvex20_over",
        showCategory = 15)
dotplot(WTvex20_under_common_GOBP,
        title = "WTvex20_under",
        showCategory = 15)
dotplot(ex9vex20_over_common_GOBP,
        title = "ex9vex20_over",
        showCategory = 15)
dotplot(ex9vex20_under_common_GOBP,
        title = "ex9vex20_under",
        showCategory = 15)
#MF
WTvex9_over_common_GOMF <- enrichGO(WTvex9_over_clean,
                                    ont = "MF",
                                    OrgDb="org.Hs.eg.db",
                                    pvalueCutoff = 0.05,
                                    keyType = "ENSEMBL",
                                    pAdjustMethod = "BH")
WTvex9_under_common_GOMF <- enrichGO(WTvex9_under_clean,
                                     ont = "MF",
                                     OrgDb="org.Hs.eg.db",
                                     pvalueCutoff = 0.05,
                                     keyType = "ENSEMBL",
                                     pAdjustMethod = "BH")
WTvex20_over_common_GOMF <- enrichGO(WTvex20_over_clean,
                                     ont = "MF",
                                     OrgDb="org.Hs.eg.db",
                                     pvalueCutoff = 0.05,
                                     keyType = "ENSEMBL",
                                     pAdjustMethod = "BH")
WTvex20_under_common_GOMF <- enrichGO(WTvex20_under_clean,
                                      ont = "MF",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENSEMBL",
                                      pAdjustMethod = "BH")
ex9vex20_over_common_GOMF <- enrichGO(ex9vex20_over_clean,
                                      ont = "MF",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENSEMBL",
                                      pAdjustMethod = "BH")
ex9vex20_under_common_GOMF <- enrichGO(ex9vex20_under_clean,
                                       ont = "MF",
                                       OrgDb="org.Hs.eg.db",
                                       pvalueCutoff = 0.05,
                                       keyType = "ENSEMBL",
                                       pAdjustMethod = "BH")
#MFplots
dotplot(WTvex9_over_common_GOMF,
        title = "WTvex9_overMF",
        showCategory = 15)
dotplot(WTvex9_under_common_GOMF,
        title = "WTvex9_underMF",
        showCategory = 15)
dotplot(WTvex20_over_common_GOMF,
        title = "WTvex20_overMF",
        showCategory = 15)
dotplot(WTvex20_under_common_GOMF,
        title = "WTvex20_underMF",
        showCategory = 15)
dotplot(ex9vex20_over_common_GOMF,
        title = "ex9vex20_overMF",
        showCategory = 15)
dotplot(ex9vex20_under_common_GOMF,
        title = "ex9vex20_underMF",
        showCategory = 15)

#combined over/under
#combine the dfs
WTvex9_combined <- c(WTvex9_under_clean, WTvex9_over_clean)
WTvex20_combined <- c(WTvex20_under_clean, WTvex20_over_clean)
ex9vex20_combined <- c(ex9vex20_under_clean, ex9vex20_over_clean)
#GOBP
WTvex9_combined_GOBP <- enrichGO(WTvex9_combined,
                                    ont = "BP",
                                    OrgDb="org.Hs.eg.db",
                                    pvalueCutoff = 0.05,
                                    keyType = "ENSEMBL",
                                    pAdjustMethod = "BH")
WTvex20_combined_GOBP <- enrichGO(WTvex20_combined,
                                 ont = "BP",
                                 OrgDb="org.Hs.eg.db",
                                 pvalueCutoff = 0.05,
                                 keyType = "ENSEMBL",
                                 pAdjustMethod = "BH")
ex9vex20_combined_GOBP <- enrichGO(ex9vex20_combined,
                                 ont = "BP",
                                 OrgDb="org.Hs.eg.db",
                                 pvalueCutoff = 0.05,
                                 keyType = "ENSEMBL",
                                 pAdjustMethod = "BH")
dotplot(WTvex9_combined_GOBP,
        title = "WTvex9_combined_GOBP",
        showCategory = 15)
dotplot(WTvex20_combined_GOBP,
        title = "WTvex20_combined_GOBP",
        showCategory = 15)
dotplot(ex9vex20_combined_GOBP,
        title = "ex9vex20_combined_GOBP",
        showCategory = 15)

#create over/under superlist without cancer cells
WTvex9_over_nocan <- intersect(MCF10A_WTvex9_over$rowname,
                                HTert_WTvex9_over$rowname)
WTvex20_over_nocan <- intersect(MCF10A_WTvex20_over$rowname,
                                 HTert_WTvex20_over$rowname)
ex9vex20_over_nocan <- intersect(MCF10A_ex9vex20_over$rowname,
                                 HTert_ex9vex20_over$rowname)
WTvex9_under_nocan <- intersect(MCF10A_WTvex9_under$rowname,
                               HTert_WTvex9_under$rowname)
WTvex20_under_nocan <- intersect(MCF10A_WTvex20_under$rowname,
                                HTert_WTvex20_under$rowname)
ex9vex20_under_nocan <- intersect(MCF10A_ex9vex20_under$rowname,
                                 HTert_ex9vex20_under$rowname)
###much more genes!
#remove transcript IDs
WTvex9_over_clean_nocan <- tools::file_path_sans_ext(
  c(WTvex9_over_nocan))
WTvex9_under_clean_nocan <- tools::file_path_sans_ext(
  c(WTvex9_under_nocan))
WTvex20_over_clean_nocan <- tools::file_path_sans_ext(
  c(WTvex20_over_nocan))
WTvex20_under_clean_nocan <- tools::file_path_sans_ext(
  c(WTvex20_under_nocan))
ex9vex20_over_clean_nocan <- tools::file_path_sans_ext(
  c(ex9vex20_over_nocan))
ex9vex20_under_clean_nocan <- tools::file_path_sans_ext(
  c(ex9vex20_under_nocan))

#venn diagram for no cancer
venn_data <- list("WTvex9_over" = c(WTvex9_over_clean_nocan),
                  "WTvex9_under" = c (WTvex9_under_clean_nocan),
                  "WTvex20_over" = c(WTvex20_over_clean_nocan),
                  "WTvex20_under" = c(WTvex20_under_clean_nocan),
                  "ex9vex20_over" = c(ex9vex20_over_clean_nocan),
                  "ex9vex20_under" = c(ex9vex20_under_clean_nocan))
p <- ggvenn(venn_data, c("WTvex9_under", "WTvex20_under",
                         "ex9vex20_over", "ex9vex20_under"),
            stroke_size = 0.5, set_name_size = 2)
p

#GOBP no cancer
WTvex9_over_common_GOBP_nocan <- enrichGO(WTvex9_over_clean_nocan,
                                    ont = "BP",
                                    OrgDb="org.Hs.eg.db",
                                    pvalueCutoff = 0.05,
                                    keyType = "ENSEMBL",
                                    pAdjustMethod = "BH")
WTvex9_under_common_GOBP_nocan <- enrichGO(WTvex9_under_clean_nocan,
                                     ont = "BP",
                                     OrgDb="org.Hs.eg.db",
                                     pvalueCutoff = 0.05,
                                     keyType = "ENSEMBL",
                                     pAdjustMethod = "BH")
WTvex20_over_common_GOBP_nocan <- enrichGO(WTvex20_over_clean_nocan,
                                     ont = "BP",
                                     OrgDb="org.Hs.eg.db",
                                     pvalueCutoff = 0.05,
                                     keyType = "ENSEMBL",
                                     pAdjustMethod = "BH")

WTvex20_under_common_GOBP_nocan <- enrichGO(WTvex20_under_clean_nocan,
                                      ont = "BP",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENSEMBL",
                                      pAdjustMethod = "BH")
ex9vex20_over_common_GOBP_nocan <- enrichGO(ex9vex20_over_clean_nocan,
                                      ont = "BP",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENSEMBL",
                                      pAdjustMethod = "BH")
ex9vex20_under_common_GOBP_nocan <- enrichGO(ex9vex20_under_clean_nocan,
                                       ont = "BP",
                                       OrgDb="org.Hs.eg.db",
                                       pvalueCutoff = 0.05,
                                       keyType = "ENSEMBL",
                                       pAdjustMethod = "BH")
dotplot(WTvex9_over_common_GOBP_nocan,
        title = "WTvex9_over",
        showCategory = 15)
dotplot(WTvex9_under_common_GOBP_nocan,
        title = "WTvex9_under",
        showCategory = 15)
dotplot(WTvex20_over_common_GOBP_nocan,
        title = "WTvex20_over",
        showCategory = 15)
dotplot(WTvex20_under_common_GOBP_nocan,
        title = "WTvex20_under",
        showCategory = 15)
dotplot(ex9vex20_over_common_GOBP_nocan,
        title = "ex9vex20_over",
        showCategory = 15)
dotplot(ex9vex20_under_common_GOBP_nocan,
        title = "ex9vex20_under",
        showCategory = 15)

#combine nocan under/over
WTvex9_combined_nocan <- c(WTvex9_under_clean_nocan, WTvex9_over_clean_nocan)
WTvex20_combined_nocan <- c(WTvex20_under_clean_nocan, WTvex20_over_clean_nocan)
ex9vex20_combined_nocan <- c(ex9vex20_under_clean_nocan, ex9vex20_over_clean_nocan)

#combined nocan GOBP
WTvex9_combined_GOBP_nocan <- enrichGO(WTvex9_combined_nocan,
                                 ont = "BP",
                                 OrgDb="org.Hs.eg.db",
                                 pvalueCutoff = 0.05,
                                 keyType = "ENSEMBL",
                                 pAdjustMethod = "BH")
WTvex20_combined_GOBP_nocan <- enrichGO(WTvex20_combined_nocan,
                                  ont = "BP",
                                  OrgDb="org.Hs.eg.db",
                                  pvalueCutoff = 0.05,
                                  keyType = "ENSEMBL",
                                  pAdjustMethod = "BH")
ex9vex20_combined_GOBP_nocan <- enrichGO(ex9vex20_combined_nocan,
                                   ont = "BP",
                                   OrgDb="org.Hs.eg.db",
                                   pvalueCutoff = 0.05,
                                   keyType = "ENSEMBL",
                                   pAdjustMethod = "BH")
dotplot(WTvex9_combined_GOBP_nocan,
        title = "WTvex9_combined_GOBP",
        showCategory = 15)
dotplot(WTvex20_combined_GOBP_nocan,
        title = "WTvex20_combined_GOBP",
        showCategory = 15)
dotplot(ex9vex20_combined_GOBP_nocan,
        title = "ex9vex20_combined_GOBP",
        showCategory = 15)
