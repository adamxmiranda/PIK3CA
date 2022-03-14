###
#Clustering on each lineage (and cancer cells) separately
#Attempt to identify genes/pathways that define mutation
#state
##Setup
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

###Cancer cell lines
#create cancer cell line counts matrix
MCF_7_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_MCF.txt", skip = 1,
                                    col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3"),
                                    row.names = "geneID")  %>%
  dplyr::select("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3")
MCF_7_counts <- MCF_7_counts %>%
  select(MCF7_corrected_rep1, MCF7_corrected_rep2, MCF7_corrected_rep3,
         MCF7_parental_rep1, MCF7_parental_rep2, MCF7_parental_rep3)
T47D_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_T47D.txt", skip = 1,
                          col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "T47D_rep1", "T47D_rep2", "T47D_rep3"),
                          row.names = "geneID") %>%
  dplyr::select("T47D_rep1", "T47D_rep2", "T47D_rep3")
cancer_counts <- cbind(MCF_7_counts, T47D_counts)
#clean data for DEseq
cancer_samples <- c("MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                    "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                    "T47D_rep1", "T47D_rep2", "T47D_rep3")
cancer_genotype <- c("WT", "WT", "WT",
                    "ex9", "ex9", "ex9",
                    "ex20", "ex20", "ex20")
cancer_genotype <- as.factor(cancer_genotype)
cancer_coldata <- data.frame(cancer_genotype)
row.names(cancer_coldata) <- cancer_samples
all(rownames(cancer_coldata) == colnames(cancer_counts))
#DEseq
cancer_dds <- DESeqDataSetFromMatrix(countData= cancer_counts,
                                    colData= cancer_coldata,
                                    design= ~ cancer_genotype)
#Filter lowcount genes and generate normalized counts
###filter = at least 20 counts in over half of the samples "5/9"
cancer_filt <- rowSums( counts(cancer_dds, normalized=FALSE) >= 20 ) >= 5

cancer_dds_f <- cancer_dds[cancer_filt,]
cancer_dds_f <- DESeq(cancer_dds_f)
#Create normalized counts file from filtered DEseq object
cancer_normalized_counts_f <- counts(cancer_dds_f, normalized=T)
View(cancer_normalized_counts_f)
cancer_normalized_counts_f_nonames <- cbind(rownames(cancer_normalized_counts_f), data.frame(cancer_normalized_counts_f, row.names=NULL))
#calc variance stabilizing transformation and regularized log to make PCA and heatmap
cancer_vst <- vst(cancer_dds_f, blind = FALSE)
head(assay(cancer_vst), 3)
cancer_rlog <- rlog(cancer_dds_f, blind = FALSE)
head(assay(cancer_rlog), 3)
#create heatmap based on sample distances of VST
cancerDists <- dist(t(assay(cancer_vst)))
#extracts the rlog matrix from the object
cancer_sampleDistMatrix <- as.matrix(cancerDists)
rownames(cancer_sampleDistMatrix) <- NULL
colnames(cancer_sampleDistMatrix) <- paste(cancer_samples)
pheatmap(cancer_sampleDistMatrix,
         clustering_distance_rows = cancerDists,
         clustering_distance_cols = cancerDists)
#PCA
plotPCA(cancer_rlog, intgroup = "cancer_genotype",)
#K-means clustering
#using the normalized counts
kmax <- 15
wss <-sapply(1:kmax,
             function(k){kmeans(cancer_normalized_counts_f, k,
                                nstart = 50, iter.max= 15) $tot.withinss})
wss
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
###using elbow method, maybe 2 or three clusters /\(`_`)/\ *big shrug*
#let's try both 2 and 3 clusters
#create heatmap with clustering
p <- pheatmap(cancer_normalized_counts_f, kmeans_k = 2)
p <- pheatmap(cancer_normalized_counts_f, kmeans_k = 3)
#honestly not sure what to make of this


###I guess I'll try the other lineages
#here we go
#MCF10A
MCF10A_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_MCF.txt", skip = 1,
                            col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3"),
                            row.names = "geneID")  %>%
  dplyr::select("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3")
MCF10A_H1047R_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_MCF10A_H1047R.txt", skip = 1,
                                   col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3"),
                                   row.names = "geneID")  %>%
  dplyr::select("MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
MCF10A_counts <- cbind(MCF10A_counts, MCF10A_H1047R_counts)
MCF10A_samples <- c("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                    "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                    "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
MCF10A_genotype <- c("ex9", "ex9", "ex9",
                     "WT", "WT", "WT",
                     "ex20", "ex20", "ex20")
MCF10A_genotype <- as.factor(MCF10A_genotype)
MCF10A_coldata <- data.frame(MCF10A_genotype)
row.names(MCF10A_coldata) <- MCF10A_samples
all(rownames(MCF10A_coldata) == colnames(MCF10A_counts))
#DEseq
MCF10A_dds <- DESeqDataSetFromMatrix(countData= MCF10A_counts,
                                     colData= MCF10A_coldata,
                                     design= ~ MCF10A_genotype)
#Filter lowcount genes and generate normalized counts
###filter = at least 20 counts in over half of the samples "5/9"
MCF10A_filt <- rowSums( counts(MCF10A_dds, normalized=FALSE) >= 20 ) >= 5
MCF10A_dds_f <- MCF10A_dds[MCF10A_filt,]
MCF10A_dds_f <- DESeq(MCF10A_dds_f)
#Create normalized counts file from filtered DEseq object
MCF10A_normalized_counts_f <- counts(MCF10A_dds_f, normalized=T)
View(MCF10A_normalized_counts_f)
MCF10A_normalized_counts_f_nonames <- cbind(rownames(MCF10A_normalized_counts_f), data.frame(MCF10A_normalized_counts_f, row.names=NULL))
#calc variance stabilizing transformation and regularized log to make PCA and heatmap
MCF10A_vst <- vst(MCF10A_dds_f, blind = FALSE)
head(assay(MCF10A_vst), 3)
MCF10A_rlog <- rlog(MCF10A_dds_f, blind = FALSE)
head(assay(MCF10A_rlog), 3)
#create heatmap based on sample distances of VST
MCF10ADists <- dist(t(assay(MCF10A_vst)))
#extracts the rlog matrix from the object
MCF10A_sampleDistMatrix <- as.matrix(MCF10ADists)
rownames(MCF10A_sampleDistMatrix) <- NULL
colnames(MCF10A_sampleDistMatrix) <- paste(MCF10A_samples)
pheatmap(MCF10A_sampleDistMatrix,
         clustering_distance_rows = MCF10ADists,
         clustering_distance_cols = MCF10ADists)
#PCA
plotPCA(MCF10A_rlog, intgroup = "MCF10A_genotype",)
#K-means clustering
#using the normalized counts
kmax <- 15
wss <-sapply(1:kmax,
             function(k){kmeans(MCF10A_normalized_counts_f, k,
                                nstart = 50, iter.max= 15) $tot.withinss})
wss
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
###using elbow method, maybe 3 or 4 clusters /\(`_`)/\ *big shrug*
#let's try both 3 and 4 clusters
#create heatmap with clustering
p <- pheatmap(MCF10A_normalized_counts_f, kmeans_k = 3)
p <- pheatmap(MCF10A_normalized_counts_f, kmeans_k = 4)
#still not sure what the hell to do with this



###I guess I'll try the last one
#here we go
#HTert
HTert_WT_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_H-Tert_WT.txt", skip = 1,
                              col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3"),
                              row.names = "geneID") %>%
  dplyr::select("HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3")
HTert_ex9_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_H-Tert_ex9.txt", skip = 1,
                               col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3"),
                               row.names = "geneID") %>%
  dplyr::select("HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3")
HTert_ex20_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_H-Tert_ex20.txt", skip = 1,
                                col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3"),
                                row.names = "geneID") %>%
  dplyr::select("HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3")
HTert_counts <- cbind(HTert_WT_counts, HTert_ex9_counts, HTert_ex20_counts)
HTert_samples <- c("HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                   "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                   "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3")
HTert_genotype <- c("WT", "WT", "WT",
                    "ex9", "ex9", "ex9",
                    "ex20", "ex20", "ex20")
HTert_genotype <- as.factor(HTert_genotype)
HTert_coldata <- data.frame(HTert_genotype)
row.names(HTert_coldata) <- HTert_samples
all(rownames(HTert_coldata) == colnames(HTert_counts))
#DEseq
HTert_dds <- DESeqDataSetFromMatrix(countData= HTert_counts,
                                     colData= HTert_coldata,
                                     design= ~ HTert_genotype)
#Filter lowcount genes and generate normalized counts
###filter = at least 20 counts in over half of the samples "5/9"
HTert_filt <- rowSums( counts(HTert_dds, normalized=FALSE) >= 20 ) >= 5
HTert_dds_f <- HTert_dds[HTert_filt,]
HTert_dds_f <- DESeq(HTert_dds_f)
#Create normalized counts file from filtered DEseq object
HTert_normalized_counts_f <- counts(HTert_dds_f, normalized=T)
View(HTert_normalized_counts_f)
HTert_normalized_counts_f_nonames <- cbind(rownames(HTert_normalized_counts_f),
                                           data.frame(HTert_normalized_counts_f,
                                                      row.names=NULL))
#calc variance stabilizing transformation and regularized log to make PCA and heatmap
HTert_vst <- vst(HTert_dds_f, blind = FALSE)
head(assay(HTert_vst), 3)
HTert_rlog <- rlog(HTert_dds_f, blind = FALSE)
head(assay(HTert_rlog), 3)
#create heatmap based on sample distances of VST
HTertDists <- dist(t(assay(HTert_vst)))
#extracts the rlog matrix from the object
HTert_sampleDistMatrix <- as.matrix(HTertDists)
rownames(HTert_sampleDistMatrix) <- NULL
colnames(HTert_sampleDistMatrix) <- paste(HTert_samples)
pheatmap(HTert_sampleDistMatrix,
         clustering_distance_rows = HTertDists,
         clustering_distance_cols = HTertDists)
#PCA
plotPCA(HTert_rlog, intgroup = "HTert_genotype",)
#K-means clustering
#using the normalized counts
kmax <- 15
wss <-sapply(1:kmax,
             function(k){kmeans(HTert_normalized_counts_f, k,
                                nstart = 50, iter.max= 15) $tot.withinss})
wss
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
###using elbow method, maybe 4 or 5 clusters /\(`_`)/\ *big shrug*
#let's try both 4 and 5 clusters
#create heatmap with clustering
p <- pheatmap(HTert_normalized_counts_f, kmeans_k = 4)
p <- pheatmap(HTert_normalized_counts_f, kmeans_k = 5)
#Ok WTF, there has to be more than like 10 genes in these clusters

###I guess we can try to cluster of diff expressed genes
#calculate diff expressed genes
#do all three comparisons in each lineage
#Cancer
cancer_WTvex9_res <- results(cancer_dds_f, contrast = c("cancer_genotype", "WT", "ex9"))
cancer_WTvex20_res <- results(cancer_dds_f, contrast = c("cancer_genotype", "WT", "ex20"))
cancer_ex9vex20_res <- results(cancer_dds_f, contrast = c("cancer_genotype", "ex9", "ex20"))
#variable name, numerator, denominator
#MCF10A
MCF10A_WTvex9_res <- results(MCF10A_dds_f, contrast = c("MCF10A_genotype", "WT", "ex9"))
MCF10A_WTvex20_res <- results(MCF10A_dds_f, contrast = c("MCF10A_genotype", "WT", "ex20"))
MCF10A_ex9vex20_res <- results(MCF10A_dds_f, contrast = c("MCF10A_genotype", "ex9", "ex20"))
#variable name, numerator, denominator
#HTert
HTert_WTvex9_res <- results(HTert_dds_f, contrast = c("HTert_genotype", "WT", "ex9"))
HTert_WTvex20_res <- results(HTert_dds_f, contrast = c("HTert_genotype", "WT", "ex20"))
HTert_ex9vex20_res <- results(HTert_dds_f, contrast = c("HTert_genotype", "ex9", "ex20"))
#variable name, numerator, denominator
#select significance based on p-adjusted values
#Cancer
cancer_WTvex9_resSig_padj <- cancer_WTvex9_res[which(cancer_WTvex9_res$padj < 0.05 ), ]
#rowcounts 11662 (nrows)
cancer_WTvex20_resSig_padj <- cancer_WTvex20_res[which(cancer_WTvex20_res$padj < 0.05 ), ]
#rowcounts 17521 (nrows)
cancer_ex9vex20_resSig_padj <- cancer_ex9vex20_res[which(cancer_ex9vex20_res$padj < 0.05 ), ]
#rowcounts 17404 (nrows)
#MCF10A
MCF10A_WTvex9_resSig_padj <- MCF10A_WTvex9_res[which(MCF10A_WTvex9_res$padj < 0.05 ), ]
#rowcounts 10199 (nrows)
MCF10A_WTvex20_resSig_padj <- MCF10A_WTvex20_res[which(MCF10A_WTvex20_res$padj < 0.05 ), ]
#rowcounts 13810 (nrows)
MCF10A_ex9vex20_resSig_padj <- MCF10A_ex9vex20_res[which(MCF10A_ex9vex20_res$padj < 0.05 ), ]
#rowcounts 13387 (nrows)
#HTert
HTert_WTvex9_resSig_padj <- HTert_WTvex9_res[which(HTert_WTvex9_res$padj < 0.05 ), ]
#rowcounts 12984 (nrows)
HTert_WTvex20_resSig_padj <- HTert_WTvex20_res[which(HTert_WTvex20_res$padj < 0.05 ), ]
#rowcounts 12047 (nrows)
HTert_ex9vex20_resSig_padj <- HTert_ex9vex20_res[which(HTert_ex9vex20_res$padj < 0.05 ), ]
#rowcounts 13967 (nrows)

#make all of these dataframes
cancer_WTvex9_resSig_padj <- as.data.frame(cancer_WTvex9_resSig_padj)
cancer_WTvex20_resSig_padj <- as.data.frame(cancer_WTvex20_resSig_padj)
cancer_ex9vex20_resSig_padj <- as.data.frame(cancer_ex9vex20_resSig_padj)
MCF10A_WTvex9_resSig_padj <- as.data.frame(MCF10A_WTvex9_resSig_padj)
MCF10A_WTvex20_resSig_padj <- as.data.frame(MCF10A_WTvex20_resSig_padj)
MCF10A_ex9vex20_resSig_padj <- as.data.frame(MCF10A_ex9vex20_resSig_padj)
HTert_WTvex9_resSig_padj <- as.data.frame(HTert_WTvex9_resSig_padj)
HTert_WTvex20_resSig_padj <- as.data.frame(HTert_WTvex20_resSig_padj)
HTert_ex9vex20_resSig_padj <- as.data.frame(HTert_ex9vex20_resSig_padj)

#filter by log2FoldChange
cancer_WTvex9_resSig_padj <- cancer_WTvex9_resSig_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
#819
cancer_WTvex20_resSig_padj <- cancer_WTvex20_resSig_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
#5431
cancer_ex9vex20_resSig_padj <- cancer_ex9vex20_resSig_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
#5192
MCF10A_WTvex9_resSig_padj <- MCF10A_WTvex9_resSig_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
#753
MCF10A_WTvex20_resSig_padj <- MCF10A_WTvex20_resSig_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
#1684
MCF10A_ex9vex20_resSig_padj <- MCF10A_ex9vex20_resSig_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
#1271
HTert_WTvex9_resSig_padj <- HTert_WTvex9_resSig_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
#1851
HTert_WTvex20_resSig_padj <- HTert_WTvex20_resSig_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
#1556
HTert_ex9vex20_resSig_padj <- HTert_ex9vex20_resSig_padj %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
#2702
#rownames to column
cancer_WTvex9_resSig_padj <- rownames_to_column(cancer_WTvex9_resSig_padj)
cancer_WTvex20_resSig_padj <- rownames_to_column(cancer_WTvex20_resSig_padj)
cancer_ex9vex20_resSig_padj <- rownames_to_column(cancer_ex9vex20_resSig_padj)
MCF10A_WTvex9_resSig_padj <- rownames_to_column(MCF10A_WTvex9_resSig_padj)
MCF10A_WTvex20_resSig_padj <- rownames_to_column(MCF10A_WTvex20_resSig_padj)
MCF10A_ex9vex20_resSig_padj <- rownames_to_column(MCF10A_ex9vex20_resSig_padj)
HTert_WTvex9_resSig_padj <- rownames_to_column(HTert_WTvex9_resSig_padj)
HTert_WTvex20_resSig_padj <- rownames_to_column(HTert_WTvex20_resSig_padj)
HTert_ex9vex20_resSig_padj <- rownames_to_column(HTert_ex9vex20_resSig_padj)

#create lists of significant genes in each lineage across all 3 comparisons
cancer_diff_genes <- rbind(cancer_WTvex9_resSig_padj, cancer_WTvex20_resSig_padj, cancer_ex9vex20_resSig_padj)
cancer_diff_genes <- cancer_diff_genes %>%
  distinct(rowname, .keep_all = TRUE)
#6503
MCF10A_diff_genes <- rbind(MCF10A_WTvex9_resSig_padj, MCF10A_WTvex20_resSig_padj, MCF10A_ex9vex20_resSig_padj)
MCF10A_diff_genes <- MCF10A_diff_genes %>%
  distinct(rowname, .keep_all = TRUE)
#2573
HTert_diff_genes <- rbind(HTert_WTvex9_resSig_padj, HTert_WTvex20_resSig_padj, HTert_ex9vex20_resSig_padj)
HTert_diff_genes <- HTert_diff_genes %>%
  distinct(rowname, .keep_all = TRUE)
#3685

#Subset the normalized gene counts for each lineage with these significant genes
#cancer
cancer_norm_sigset <- subset(cancer_normalized_counts_f_nonames, cancer_normalized_counts_f_nonames$`rownames(cancer_normalized_counts_f)` %in%
                               cancer_diff_genes$rowname)
cancer_norm_sigset <- remove_rownames(cancer_norm_sigset)
cancer_norm_sigset <- as.data.frame(cancer_norm_sigset)
cancer_norm_sigset <- column_to_rownames(cancer_norm_sigset,
                                         var = "rownames(cancer_normalized_counts_f)")
#MCF10A
MCF10A_norm_sigset <- subset(MCF10A_normalized_counts_f_nonames, MCF10A_normalized_counts_f_nonames$`rownames(MCF10A_normalized_counts_f)` %in%
                               MCF10A_diff_genes$rowname)
MCF10A_norm_sigset <- remove_rownames(MCF10A_norm_sigset)
MCF10A_norm_sigset <- as.data.frame(MCF10A_norm_sigset)
MCF10A_norm_sigset <- column_to_rownames(MCF10A_norm_sigset,
                                         var = "rownames(MCF10A_normalized_counts_f)")
#HTert
HTert_norm_sigset <- subset(HTert_normalized_counts_f_nonames, HTert_normalized_counts_f_nonames$`rownames(HTert_normalized_counts_f)` %in%
                              HTert_diff_genes$rowname)
HTert_norm_sigset <- remove_rownames(HTert_norm_sigset)
HTert_norm_sigset <- as.data.frame(HTert_norm_sigset)
HTert_norm_sigset <- column_to_rownames(HTert_norm_sigset,
                                         var = "rownames(HTert_normalized_counts_f)")

#k-means with the subsets
#Cancer
kmax <- 15
wss <-sapply(1:kmax,
             function(k){kmeans(cancer_norm_sigset, k,
                                nstart = 50, iter.max= 15) $tot.withinss})
wss
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
p <- pheatmap(cancer_norm_sigset)
#MCF10A
wss <-sapply(1:kmax,
             function(k){kmeans(MCF10A_norm_sigset, k,
                                nstart = 50, iter.max= 15) $tot.withinss})
wss
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
p <- pheatmap(MCF10A_norm_sigset)
#HTert
wss <-sapply(1:kmax,
             function(k){kmeans(HTert_norm_sigset, k,
                                nstart = 50, iter.max= 15) $tot.withinss})
wss
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
p <- pheatmap(HTert_norm_sigset, kmeans_k = 3)
