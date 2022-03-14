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

#selected 17635 genes with padj < 0.05 
MCF10A_ex20vparental_genes_padj <- as.data.frame(read.csv("MCF10A_H1047Rvparental.csv"))
MCF10A_ex20vparental_genes_padj_entrez <- as.character(MCF10A_ex20vparental_genes_padj$entrez)

#selected 16924 genes with padj < 0.05 
MCF10A_ex20vex9_genes_padj <- as.data.frame(read.csv("MCF10A_H1047RvE545K.csv"))
MCF10A_ex20vex9_genes_padj_entrez <- as.character(MCF10A_ex20vex9_genes_padj$entrez)

#selected 10256 genes with padj < 0.05 
MCF10A_ex9vparental_genes_padj <- as.data.frame(read.csv("MCF10_E545Kvparental.csv"))
MCF10A_ex9vparental_genes_padj_entrez <- as.character(MCF10A_ex9vparental_genes_padj$entrez)

#filter by l2fc
#selected 5650 genes with log2foldchange outside of (+/-)1.5 and padj < 0.05 
MCF10A_ex20vparental_genes_padj_log2fc_1.5 <- MCF10A_ex20vparental_genes_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
MCF10A_ex20vparental_genes_padj_log2fc_1.5_entrez <- as.character(MCF10A_ex20vparental_genes_padj_log2fc_1.5$entrez)
MCF10A_ex20vparental_genes_padj_log2fc_1.5_gene_list_df<- as.data.frame(MCF10A_ex20vparental_genes_padj_log2fc_1.5$entrez)
write.csv(MCF10A_ex20vparental_genes_padj_log2fc_1.5, file = "MCF10A_ex20vparental_Comparison_Sig_entrez_list")

#selected 4998 genes with log2foldchange outside of (+/-)1.5 and padj < 0.05
MCF10A_ex20vex9_genes_padj_log2fc_1.5 <- MCF10A_ex20vex9_genes_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
MCF10A_ex20vex9_genes_padj_log2fc_1.5_entrez <- as.character(MCF10A_ex20vex9_genes_padj_log2fc_1.5$entrez)
MCF10A_ex20vex9_genes_padj_log2fc_1.5_entrez_df <- as.data.frame(MCF10A_ex20vex9_genes_padj_log2fc_1.5$entrez)
write.csv(MCF10A_ex20vex9_genes_padj_log2fc_1.5, file = "MCF10_ex20vex9_Comparison_Sig_entrez_list")

#selected 1276 genes with log2foldchange outside of (+/-)1.5 and padj < 0.05
MCF10A_ex9vparental_genes_padj_log2fc_1.5 <- MCF10A_ex9vparental_genes_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
MCF10A_ex9vparental_genes_padj_log2fc_1.5_entrez <- as.character(MCF10A_ex9vparental_genes_padj_log2fc_1.5$entrez)
MCF10A_ex9vparental_genes_padj_log2fc_1.5_entrez_df <- as.data.frame(MCF10A_ex9vparental_genes_padj_log2fc_1.5$entrez)
write.csv(MCF10A_ex9vparental_genes_padj_log2fc_1.5, file = "MCF10_ex9vparental_Comparison_Sig_entrez_list")

######GO Biological Process######
MCF10A_ex20vparental_GOBP <- enrichGO(MCF10A_ex20vparental_genes_padj_log2fc_1.5_entrez,
                       ont = "BP",
                       OrgDb="org.Hs.eg.db",
                       pvalueCutoff = 0.05,
                       keyType = "ENTREZID",
                       pAdjustMethod = "BH")
dotplot(MCF10A_ex20vparental_GOBP,
        title = "MCF10A H1047R v parental Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

MCF10A_ex20vex9_GOBP <- enrichGO(MCF10A_ex20vex9_genes_padj_log2fc_1.5_entrez,
                                      ont = "BP",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENTREZID",
                                      pAdjustMethod = "BH")
dotplot(MCF10A_ex20vex9_GOBP,
        title = "MCF10A H1047R v E545K Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

MCF10A_ex9vparental_GOBP <- enrichGO(MCF10A_ex9vparental_genes_padj_log2fc_1.5_entrez,
                                      ont = "BP",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENTREZID",
                                      pAdjustMethod = "BH")
dotplot(MCF10A_ex9vparental_GOBP,
        title = "MCF10A E545K v parental Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

######GO Molecular Function######
MCF10A_ex20vparental_GOMF <- enrichGO(MCF10A_ex20vparental_genes_padj_log2fc_1.5_entrez,
                                      ont = "MF",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENTREZID",
                                      pAdjustMethod = "BH")
dotplot(MCF10A_ex20vparental_GOMF,
        title = "MCF10A H1047R v parental Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 15)

MCF10A_ex20vex9_GOMF <- enrichGO(MCF10A_ex20vex9_genes_padj_log2fc_1.5_entrez,
                                 ont = "MF",
                                 OrgDb="org.Hs.eg.db",
                                 pvalueCutoff = 0.05,
                                 keyType = "ENTREZID",
                                 pAdjustMethod = "BH")
dotplot(MCF10A_ex20vex9_GOMF,
        title = "MCF10A H1047R v E545K Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 20)

MCF10A_ex9vparental_GOMF <- enrichGO(MCF10A_ex9vparental_genes_padj_log2fc_1.5_entrez,
                                     ont = "MF",
                                     OrgDb="org.Hs.eg.db",
                                     pvalueCutoff = 0.05,
                                     keyType = "ENTREZID",
                                     pAdjustMethod = "BH")
dotplot(MCF10A_ex9vparental_GOMF,
        title = "MCF10A E545K v parental Padj < 0.05 & l2fc <> +/-1.5 Ontology",
        showCategory = 20)

MCF10A_ex20vparental_channel_list<- dplyr::filter(as.data.frame(MCF10A_ex20vparental_GOMF), Description == "channel activity")
MCF10A_ex20vex9_channel_list<- dplyr::filter(as.data.frame(MCF10A_ex20vex9_GOMF), Description == "channel activity")
print(MCF10A_ex20vparental_channel_list)
print(MCF10A_ex20vex9_channel_list)

MCF10A_ex20vparental_transcription_activator_RNAP2_list<- dplyr::filter(as.data.frame(MCF10A_ex20vparental_GOMF), Description == "DNA-binding transcription activator activity, RNA polymerase II-specific")
print(MCF10A_ex20vparental_transcription_activator_RNAP2_list)
MCF10A_ex20vex9_transcription_activator_RNAP2_list<- dplyr::filter(as.data.frame(MCF10A_ex20vex9_GOMF), Description == "DNA-binding transcription activator activity, RNA polymerase II-specific")
print(MCF10A_ex20vex9_transcription_activator_RNAP2_list$geneID)
TFs_shared<-intersect(MCF10A_ex20vparental_transcription_activator_RNAP2_list, MCF10A_ex20vex9_transcription_activator_RNAP2_list)
###Enrich Kegg###
MCF10A_ex20vparental_KEGG <- enrichKEGG(MCF10A_ex20vparental_genes_padj_log2fc_1.5_entrez,
                                         organism = "hsa",
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH")

dotplot(MCF10A_ex20vparental_KEGG,
        showCategory = 20,
        title = "L2fc(+/-1.5) KEGG")

MCF10A_ex20vparental_FAK_list <- dplyr::filter(as.data.frame(MCF10A_ex20vparental_KEGG), Description == "Focal adhesion")
print(MCF10A_ex20vparental_FAK_list$geneID)
MCF10A_ex20vparental_MAPK_list <- dplyr::filter(as.data.frame(MCF10A_ex20vparental_KEGG), Description == "MAPK signaling pathway")
print(MCF10A_ex20vparental_MAPK_list$geneID)
MCF10A_ex20vparental_PI3K_list <- "2335/894/7057/2264/1944/1280/1293/10161/3481/54331/2791/8515/1298/2784/2069/5228/7450/2324/2246/5649/3674/3667/596/118788/2277/5155/2321/22798/131873/3690"
MCF10A_ex20vparental_breast_list <- "948/348/29116/3949/255738/4036"


MCF10A_ex20vex9_KEGG <- enrichKEGG(MCF10A_ex20vex9_genes_padj_log2fc_1.5_entrez,
                                        organism = "hsa",
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")

dotplot(MCF10A_ex20vex9_KEGG,
        showCategory = 20,
        title = "L2fc(+/-1.5) KEGG")

MCF10A_ex9vparental_KEGG <- enrichKEGG(MCF10A_ex9vparental_genes_padj_log2fc_1.5_entrez,
                                   organism = "hsa",
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH")

dotplot(MCF10A_ex9vparental_KEGG,
        showCategory = 20,
        title = "L2fc(+/-1.5) KEGG")


######Interrogate very low expression genes from ex20vparental comparison######
#selected 5650 genes with log2foldchange outside of (-)5 and padj < 0.05 
MCF10A_ex20vparental_genes_padj_log2fc_neg5 <- MCF10A_ex20vparental_genes_padj %>%
  filter(log2FoldChange < -5)
MCF10A_ex20vparental_genes_padj_log2fc_neg5_entrez <- as.character(MCF10A_ex20vparental_genes_padj_log2fc_neg5$entrez)
MCF10A_ex20vparental_genes_padj_log2fc_neg5_gene_list_df<- as.data.frame(MCF10A_ex20vparental_genes_padj_log2fc_neg5$entrez)
write.csv(MCF10A_ex20vparental_genes_padj_log2fc_neg5, file = "MCF10A_ex20vparental_Comparison_Sig_neg5_entrez_list")

MCF10A_ex20vparental_neg5_GOMF <- enrichGO(MCF10A_ex20vparental_genes_padj_log2fc_neg5_entrez,
                                      ont = "MF",
                                      OrgDb="org.Hs.eg.db",
                                      pvalueCutoff = 0.05,
                                      keyType = "ENTREZID",
                                      pAdjustMethod = "BH")
dotplot(MCF10A_ex20vparental_neg5_GOMF,
        title = "MCF10A H1047R v parental Padj < 0.05 & l2fc -5 Ontology",
        showCategory = 15)

MCF10A_ex20vparental_neg5_GOBP <- enrichGO(MCF10A_ex20vparental_genes_padj_log2fc_neg5_entrez,
                                           ont = "BP",
                                           OrgDb="org.Hs.eg.db",
                                           pvalueCutoff = 0.05,
                                           keyType = "ENTREZID",
                                           pAdjustMethod = "BH")
dotplot(MCF10A_ex20vparental_neg5_GOBP,
        title = "MCF10A H1047R v parental Padj < 0.05 & l2fc -5 Ontology",
        showCategory = 15)

MCF10A_ex20vparental_transcription_activator_RNAP2_list_neg5<- dplyr::filter(as.data.frame(MCF10A_ex20vparental_neg5_GOMF), Description == "DNA-binding transcription activator activity, RNA polymerase II-specific")
print(MCF10A_ex20vparental_transcription_activator_RNAP2_list_neg5)
######Interrogate very upregulated genes from ex20vparental comparison######
#selected 1794 genes with log2foldchange outside of +1.5 and padj < 0.05 
MCF10A_ex20vparental_genes_padj_log2fc_upreg <- MCF10A_ex20vparental_genes_padj %>%
  filter(log2FoldChange > 1.5)
MCF10A_ex20vparental_genes_padj_log2fc_upreg_entrez <- as.character(MCF10A_ex20vparental_genes_padj_log2fc_upreg$entrez)
MCF10A_ex20vparental_genes_padj_log2fc_upreg_gene_list_df<- as.data.frame(MCF10A_ex20vparental_genes_padj_log2fc_upreg$entrez)
write.csv(MCF10A_ex20vparental_genes_padj_log2fc_upreg, file = "MCF10A_ex20vparental_Comparison_Sig_upreg_entrez_list")

MCF10A_ex20vparental_upreg_GOMF <- enrichGO(MCF10A_ex20vparental_genes_padj_log2fc_upreg_entrez,
                                           ont = "MF",
                                           OrgDb="org.Hs.eg.db",
                                           pvalueCutoff = 0.05,
                                           keyType = "ENTREZID",
                                           pAdjustMethod = "BH")
dotplot(MCF10A_ex20vparental_upreg_GOMF,
        title = "MCF10A H1047R v parental Padj < 0.05 & l2fc +1.5 Ontology",
        showCategory = 15)

MCF10A_ex20vparental_upreg_GOBP <- enrichGO(MCF10A_ex20vparental_genes_padj_log2fc_upreg_entrez,
                                           ont = "BP",
                                           OrgDb="org.Hs.eg.db",
                                           pvalueCutoff = 0.05,
                                           keyType = "ENTREZID",
                                           pAdjustMethod = "BH")
dotplot(MCF10A_ex20vparental_upreg_GOBP,
        title = "MCF10A H1047R v parental Padj < 0.05 & l2fc +1.5 Ontology",
        showCategory = 15)

MCF10A_ex20vparental_upreg_KEGG <- enrichKEGG(MCF10A_ex20vparental_genes_padj_log2fc_upreg_entrez,
                                         organism = "hsa",
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH")

dotplot(MCF10A_ex20vparental_upreg_KEGG,
        showCategory = 20,
        title = "L2fc(+/-1.5) KEGG")


######Interrogate very upregulated genes from ex20vex9 comparison######
#selected 1260 genes with log2foldchange outside of +1.5 and padj < 0.05 
MCF10A_ex20vex9_genes_padj_log2fc_upreg <- MCF10A_ex20vex9_genes_padj %>%
  filter(log2FoldChange > 1.5)
MCF10A_ex20vex9_genes_padj_log2fc_upreg_entrez <- as.character(MCF10A_ex20vex9_genes_padj_log2fc_upreg$entrez)
MCF10A_ex20vex9_genes_padj_log2fc_upreg_gene_list_df<- as.data.frame(MCF10A_ex20vex9_genes_padj_log2fc_upreg$entrez)
write.csv(MCF10A_ex20vex9_genes_padj_log2fc_upreg, file = "MCF10A_ex20vex9_Comparison_Sig_upreg_entrez_list")

MCF10A_ex20vex9_upreg_GOMF <- enrichGO(MCF10A_ex20vex9_genes_padj_log2fc_upreg_entrez,
                                            ont = "MF",
                                            OrgDb="org.Hs.eg.db",
                                            pvalueCutoff = 0.05,
                                            keyType = "ENTREZID",
                                            pAdjustMethod = "BH")
dotplot(MCF10A_ex20vex9_upreg_GOMF,
        title = "MCF10A H1047R v E545K Padj < 0.05 & l2fc +1.5 Ontology",
        showCategory = 15)

MCF10A_ex20vex9_upreg_GOBP <- enrichGO(MCF10A_ex20vex9_genes_padj_log2fc_upreg_entrez,
                                            ont = "BP",
                                            OrgDb="org.Hs.eg.db",
                                            pvalueCutoff = 0.05,
                                            keyType = "ENTREZID",
                                            pAdjustMethod = "BH")
dotplot(MCF10A_ex20vex9_upreg_GOBP,
        title = "MCF10A H1047R v E545K Padj < 0.05 & l2fc +1.5 Ontology",
        showCategory = 15)

MCF10A_ex20vex9_upreg_KEGG <- enrichKEGG(MCF10A_ex20vex9_genes_padj_log2fc_upreg_entrez,
                                    organism = "hsa",
                                    pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH")

dotplot(MCF10A_ex20vex9_upreg_KEGG,
        showCategory = 20,
        title = "L2fc(+/-1.5) KEGG")

###Intersection of differentially expressed genes against parental from both mutants
intersection_mutantsvparental_gene_list <- intersect(MCF10A_ex20vparental_genes_padj_log2fc_1.5_entrez, MCF10A_ex9vparental_genes_padj_log2fc_1.5_entrez)
intersection_mutantsvparental_gene_list_df <- as.data.frame(intersection_mutantsvparental_gene_list)
#443 genes intersect

mutant_intersect_GOBP <- enrichGO(intersection_mutantsvparental_gene_list,
                                       ont = "BP",
                                       OrgDb="org.Hs.eg.db",
                                       pvalueCutoff = 0.05,
                                       keyType = "ENTREZID",
                                       pAdjustMethod = "BH")
dotplot(mutant_intersect_GOBP,
        title = "Mutants v parental Padj < 0.05 & l2fc +1.5 Ontology",
        showCategory = 15)

mutant_intersect_GOMF <- enrichGO(intersection_mutantsvparental_gene_list,
                                  ont = "MF",
                                  OrgDb="org.Hs.eg.db",
                                  pvalueCutoff = 0.05,
                                  keyType = "ENTREZID",
                                  pAdjustMethod = "BH")
dotplot(mutant_intersect_GOMF,
        title = "Mutants v parental Padj < 0.05 & l2fc +1.5 Ontology",
        showCategory = 15)

mutant_intersect_KEGG <- enrichKEGG(intersection_mutantsvparental_gene_list,
                                         organism = "hsa",
                                         pvalueCutoff = 1.0,
                                         pAdjustMethod = "BH")

dotplot(mutant_intersect_KEGG,
        showCategory = 20,
        title = "L2fc(+/-1.5) KEGG")

###Subset off the bump of low expression###
MCF10A_ex20vparental_genes_padj_log2fc_nobump <- MCF10A_ex20vparental_genes_padj %>%
  filter(between(log2FoldChange, -5, -1.5) | log2FoldChange > 1.5)
MCF10A_ex20vparental_genes_padj_log2fc_nobump_entrez <- as.character(MCF10A_ex20vparental_genes_padj_log2fc_nobump$entrez)
MCF10A_ex20vparental_genes_padj_log2fc_nobump_gene_list_df<- as.data.frame(MCF10A_ex20vparental_genes_padj_log2fc_nobump$entrez)
write.csv(MCF10A_ex20vparental_genes_padj_log2fc_nobump, file = "MCF10A_ex20vparental_Comparison_Sig_nobump_entrez_list")

MCF10A_ex20vparental_nobump_GOBP <- enrichGO(MCF10A_ex20vparental_genes_padj_log2fc_nobump_entrez,
                                  ont = "BP",
                                  OrgDb="org.Hs.eg.db",
                                  pvalueCutoff = 0.05,
                                  keyType = "ENTREZID",
                                  pAdjustMethod = "BH")
dotplot(MCF10A_ex20vparental_nobump_GOBP,
        title = "Mutants v parental Padj < 0.05 & l2fc +1.5 Ontology",
        showCategory = 15)

MCF10A_ex20vparental_nobump_GOMF <- enrichGO(MCF10A_ex20vparental_genes_padj_log2fc_nobump_entrez,
                                  ont = "MF",
                                  OrgDb="org.Hs.eg.db",
                                  pvalueCutoff = 0.05,
                                  keyType = "ENTREZID",
                                  pAdjustMethod = "BH")
dotplot(MCF10A_ex20vparental_nobump_GOMF,
        title = "Mutants v parental Padj < 0.05 & l2fc +1.5 Ontology",
        showCategory = 15)

MCF10A_ex20vparental_nobump_KEGG <- enrichKEGG(MCF10A_ex20vparental_genes_padj_log2fc_nobump_entrez,
                                    organism = "hsa",
                                    pvalueCutoff = 1.0,
                                    pAdjustMethod = "BH")

dotplot(MCF10A_ex20vparental_nobump_KEGG,
        showCategory = 20,
        title = "L2fc(+/-1.5) KEGG")
