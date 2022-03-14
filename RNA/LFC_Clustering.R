#trying to cluster the data the way Kelly did
#using the LFC as the plotted values rather
#than the normalized counts
library(DESeq2)
library(ashr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
#Using the total unfiltered dds DESeq object from Total_DESeq.R, line 78-81
total_dds
#standard differential expression analysis
total_dds <- DESeq(total_dds)
resultsNames(total_dds)

# Apply likelihood ratio test (LRT) to data, this design assumes
#only PIK3CA mutation as a variant
#using total_treatment from Total_DESeq.R, line 56-64
dds_LRT <- DESeq(total_dds, test="LRT",
                 full = ~total_treatment,
                 reduced = ~1)
# Write results table of LRT test. Will contain Log2 fold change, p-values, and adjusted p-values for each gene
# Genes that have too little data to make a conclusion will be assigned a value of NA. This can be changed in the results() settings
resLRT <- results(dds_LRT)

# Log transforms the DESeq2 count table so that it is more appropriate for use in clustering or PCA
rlog_dds <- rlog(total_dds, blind=F) 

# Shrink the lfc of low count genes so they do not artificially skew results 
# These are shrunken LFC in reference to different comparisons
resAshWTvex9 <- lfcShrink(dds_LRT, coef=1, type="ashr")
resAshex9vex20 <- lfcShrink(dds_LRT, coef=2, type="ashr")
resAshWTvex20 <- lfcShrink(dds_LRT, coef=3, type="ashr")

# Build new dataframe of shrunken log fold changes
lfcShrink_df <- data.frame(rownames(dds_LRT),
                           as.numeric("0"),
                           resAshWTvex9$log2FoldChange,
                           resAshex9vex20$log2FoldChange,
                           resAshWTvex20$log2FoldChange,
                           resLRT$padj)
colnames(lfcShrink_df) <- c("ENSEMBL",
                            "lfc0",
                            "lfcWTvex9",
                            "lfcex9vex20",
                            "lfcWTvex20",
                            "padj")
# Build similar dataframe but take the absolute value of log fold change
# Also add column that is the rowmax for absolute value lfc values 
lfcShrink_df_absvalue <- data.frame(rownames(dds_LRT),
                                    as.numeric("0"),
                                    abs(resAshWTvex9$log2FoldChange),
                                    abs(resAshex9vex20$log2FoldChange),
                                    abs(resAshWTvex20$log2FoldChange),
                                    resLRT$padj)
lfcShrink_df_absvalue <- data.frame(rownames(dds_LRT),
                                    as.numeric("0"),
                                    abs(resAshWTvex9$log2FoldChange),
                                    abs(resAshex9vex20$log2FoldChange),
                                    abs(resAshWTvex20$log2FoldChange),
                                    apply(lfcShrink_df_absvalue[ ,2:5],
                                          1, max),
                                    resLRT$padj)
colnames(lfcShrink_df_absvalue) <- c("ENSEMBL",
                                     "lfc0",
                                     "lfcWTvex9",
                                     "lfcex9vex20",
                                     "lfcWTvex20",
                                     "lfcMax",
                                     "padj")
# Filter lfcShrink_df_absvalue for padj value to remove low confidence genes 
lfcShrink_df_absvalue <- filter(lfcShrink_df_absvalue, padj <= .005)
# Sort data by abs(logFC) maxima
lfcShrink_df_absvalue <- arrange(lfcShrink_df_absvalue, desc(lfcMax))
# Filter for logFC >=2 between any consecutive comparisons 
lfcShrink_df_absvalue_filt <- filter(lfcShrink_df_absvalue, lfcMax >= 2)

# Filters the log transformed data against the lfc Maxima
rld.sigtest2 <- rlog_dds[ which(lfcShrink_df_absvalue$lfcMax >= 2), ]

# Generate matrix that we can feed to pheatmap for clustering
mat2 <- assay(rld.sigtest2)
# Remove rows that have no counts
mat2.1 <- mat2[ which(rowSums2(mat2) != 0), ]
# Re-order rows of  matrix for input into pheatmap
list <- c("MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
          "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
          "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3",
          "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
          "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
          "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
          "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
          "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
          "T47D_rep1", "T47D_rep2", "T47D_rep3")
kmax <- 15
wss <-sapply(1:kmax,
             function(k){kmeans(mat2.1, k,
                                nstart = 50, iter.max= 15) $tot.withinss})
wss
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
heatmapres2.1<- pheatmap(mat2.1[,list], kmeans_k = 9)

#pull genes in cluster 9
cluster_df <- as.data.frame(heatmapres2.1$kmeans$cluster)
cluster9 <- cluster_df %>%
  subset(`heatmapres2.1$kmeans$cluster` ==9)
cluster9 <- tibble::rownames_to_column(cluster9, var = "gene")
cluster9_genes <- as.vector (cluster9$gene)
View(cluster9_genes)
cluster9_genes

#pull genes in cluster 3
cluster3 <- cluster_df %>%
  subset(`heatmapres2.1$kmeans$cluster` ==3)
cluster3 <- tibble::rownames_to_column(cluster3, var = "gene")
cluster3_genes <- as.vector (cluster3$gene)
cluster3_genes

#Try GO
#remove transcript IDs
cluster9_genes_clean <- tools::file_path_sans_ext(
  c(cluster9_genes))

library("clusterProfiler")
cluster9_GOBP <- enrichGO(cluster9_genes_clean,
                          ont = "BP",
                          OrgDb="org.Hs.eg.db",
                          pvalueCutoff = 0.1,
                          keyType = "ENSEMBL",
                          pAdjustMethod = "BH")
dotplot(cluster9_GOBP,
        title = "cluster9_GOBP",
        showCategory = 15)

#cluster 3
cluster3_genes_clean <- tools::file_path_sans_ext(
  c(cluster3_genes))

cluster3_GOBP <- enrichGO(cluster3_genes_clean,
                          ont = "BP",
                          OrgDb="org.Hs.eg.db",
                          pvalueCutoff = 0.1,
                          keyType = "ENSEMBL",
                          pAdjustMethod = "BH")
dotplot(cluster3_GOBP,
        title = "cluster3_GOBP",
        showCategory = 15)

#Hierarchical Clustering
datares_HC <- pheatmap(mat2[,list], cluster_cols = T)
