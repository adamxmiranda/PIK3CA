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

setwd('~/Desktop/RNAseq_Dual')
#read in the feature counts file
#read.delim is used to read in delimited text files, where data is organized in a data matrix with rows representing cases and columns representing variables
#select() allows you to zoom in on columns that are actually of interest to you in the feature counts file using dyplyr
MCF_dual_E545K_counts <- read.delim(file = "~/Desktop/RNAseq_dual/featureCounts_all_samples_MCF.txt", skip = 1,
                              col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3"),
                              row.names = "geneID")  %>%
  dplyr::select("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3")
MCF10A_H1047R_counts <- read.delim(file = "~/Desktop/RNAseq_dual/featureCounts_all_samples_MCF10A_H1047R.txt", skip = 1,
                              col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3"),
                              row.names = "geneID")  %>%
  dplyr::select("MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
#combine dataframes with cbind
total_counts <- cbind(MCF_dual_E545K_counts, MCF10A_H1047R_counts)
MCF10A_counts <- total_counts[c("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")]
#creating a coldata table for the DESeq matrix, which is just a table with metadata on the count table's columns; colData = column metadata table
#colData = metadata = the "columns" of the count table
#colData is a data.frame that can contain all the variables you know about your samples, such as the experimental condition, the type and date of sequencing and so on
MCF_total_samples <- c("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3", "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
MCF_total_treatment<- c("10AE545K", "10AE545K", "10AE545K", "10AParental", "10AParental", "10AParental", "7Corrected", "7Corrected", "7Corrected", "7Parental","7Parental", "7Parental", "10AH1047R", "10AH1047R", "10AH1047R")
#values need to be in factor format for DESeq
MCF_total_treatment <- as.factor(MCF_total_treatment)

#character variables passed to data.frame are converted to factor columns (i.e., CRISPRi treatment)
#row names are CRISPRi samples from colData data frame
MCF_total_coldata <- data.frame(MCF_total_treatment)
row.names(MCF_total_coldata) <- MCF_total_samples
all(rownames(MCF_total_coldata) == colnames(total_counts))

#creating DESeq data object from the matrix of counts and the metadata table (colData) using DESeqDataSetFromMatrix input function
#testing for the effect of treatment (e.g., Control CRISPRi versus CRISPRi targeting ZNF384) with design variable CRISPRi_treatment
#countData = siRNA_counts = a table with the full read counts
#to avoid model matrix is not a full rank error can only consider the treatment variable in the design
#design formula accounts for source of variation (i.e., CRISPRi_treatment)
MCF_total_dds <- DESeqDataSetFromMatrix(countData= total_counts,
                                       colData= MCF_total_coldata,
                                       design= ~ MCF_total_treatment)

###Additional Pre-filtering step for low count genes
filt <- rowSums( counts(MCF_total_dds, normalized=FALSE) >= 20 ) >= 7

MCF_total_dds_f <- MCF_total_dds[filt,]
MCF_total_dds_f <- DESeq(MCF_total_dds_f)

normalized_MCF_total_counts_f <- counts(MCF_total_dds_f, normalized=T)
View(normalized_MCF_total_counts_f)

normalized_MCF_total_counts_f_nonames <- cbind(rownames(normalized_MCF_total_counts_f), data.frame(normalized_MCF_total_counts_f, row.names=NULL))
###See if this impacts the 0.00000e+00 issue
###It does not even at a 30 count threshold
###However increasing the sample threshold significantly reduces their presence
###20, 7 seems to remove all of these instances

#pre-filtering the counts so to eliminates genes that had no reads or less than 10 reads map to them
#kick out genes that have a total read count less than 10 across all samples, which leaves 34982 genes out of the total count
#pre-filtering low gene counts reduces the memory size of the dds data object (CRISPRi_dds) and increases the speed of the transformation and testing functions within DESeq2
#MCF_total_counts_filtered <- rowSums(counts(MCF_total_dds)) > 10
#MCF_total_dds <- MCF_total_dds[MCF_total_counts_filtered,]
#rowcount 44117

###TRACKSWITCH### Using the more stringent filter
#To run the differential expression analysis, we use a single call function DESeq()
#This function prints out a message for various steps: estimating pre-size factors, dispersons etc...
#has every comparison that can be drawn
MCF_total_dds <- DESeq(MCF_total_dds_f)
View(counts(MCF_total_dds))
#to see what comparisons were made; resultsNames returns the names of the estimated effects (coefficients) of the model
resultsNames(MCF_total_dds)
#To normalize the count data, DESeq2 calculates size factors (i.e., normalization factors) for each samples using the median of ratios method, which is automatically done when performing the differential expression analysis -- DESeq()
#Normalization takes into account: sequencing depth, gene length, and RNA composition
#Median of ratios method makes the assumptions that not ALL genes are differentially expressed
#Check the size factors; each sample has a unique normalization (size) factor
sizeFactors(MCF_total_dds)
#Total number of raw counts per sample
#How do the numbers correlate with size factor?
colSums(counts(MCF_total_dds))
#Total number of normalized counts per sample
#Note: DESeq2 doesn't actually used normalized counts; need raw counts as input (i.e., featureCounts txt file)
colSums(counts(MCF_total_dds, normalized=T))
normalized_MCF_total_counts <- counts(MCF_total_dds, normalized=T)

#save dds file
save(MCF_total_dds, file = "MCF_total_dds.Rdata")

#independent filtering is performed automatically with the results function
#outputting the filtered results, which extracts a results table with log2 fold changes, p values and adjusted p values
#in a previous step, we filtered out the genes that that zero counts in all samples or less than 10 reads across
#now with results() we  remove genes with an extreme count outlier and/or a low mean normalized count
#reference level = Control_CRISPRi
#e.g., log2(Parental/Corrected CRISPRi)
#this is where the comparisons are made, list all of the comparisons that are relevant to the data
MCF7_res <- results(MCF_total_dds, contrast = c("MCF_total_treatment", "7Corrected", "7Parental"))
MCF10A_E545K_res <- results(MCF_total_dds, contrast = c("MCF_total_treatment", "10AE545K", "10AParental"))
MCF_E545K_res <- results(MCF_total_dds, contrast = c("MCF_total_treatment", "10AE545K", "7Parental"))
MCF_WT_res <- results(MCF_total_dds, contrast = c("MCF_total_treatment", "10AParental", "7Corrected"))
#New comparisons with the H1047R data
MCF10A_H1047R_res <- results(MCF_total_dds, contrast = c("MCF_total_treatment", "10AH1047R", "10AParental"))
MCF10A_Mut_res <- results(MCF_total_dds, contrast = c("MCF_total_treatment", "10AH1047R", "10AE545K"))
#variable name, numerator, denominator

MCF7_res_df <- as.data.frame(MCF7_res)
MCF10A_E545K_res_df <- as.data.frame(MCF10A_E545K_res)
MCF_E545K_res_df <- as.data.frame(MCF_E545K_res)
MCF_WT_res_df <- as.data.frame(MCF_WT_res)
MCF10A_H1047R_res_df <- as.data.frame(MCF10A_H1047R_res)
MCF10A_Mut_res_df <- as.data.frame(MCF10A_Mut_res)
#data.frame allows you to use CRISPRi results for other R packages

diff_7 <- rownames_to_column(MCF7_res_df)
diff_10A_E545K <- rownames_to_column(MCF10A_E545K_res_df)
diff_E545K_Mutant <- rownames_to_column(MCF_E545K_res_df)
diff_WT <- rownames_to_column(MCF_WT_res_df)
diff_10A_H1047R <- rownames_to_column(MCF10A_H1047R_res_df)
diff_10A_Mutants <- rownames_to_column(MCF10A_Mut_res_df)
#row names become a column

mcols(MCF7_res, use.names=TRUE)
mcols(MCF10A_E545K_res, use.names=TRUE)
mcols(MCF_E545K_res, use.names=TRUE)
mcols(MCF_WT_res, use.names=TRUE)
mcols(MCF10A_H1047R_res, use.names=TRUE)
mcols(MCF10A_Mut_res, use.names=TRUE)
#carries metadata with information on the meaning of the columns
#dif_2 <- diff %>% dplyr::filter(str_detect(diff$rowname, pattern = "ENSG00000126746")) --> ENSEMBLE code for ZNF384

MCF7_resSig_padj <- MCF7_res[which(MCF7_res$padj < 0.05 ), ]
#rowcounts 12682 (nrows)
MCF10A_E545K_resSig_padj <- MCF10A_E545K_res[which(MCF10A_E545K_res$padj < 0.05 ), ]
#rowcounts 10256 (nrows)
MCF_E545K_resSig_padj <- MCF_E545K_res[which(MCF_E545K_res$padj < 0.05 ), ]
#rowcounts 22015 (nrows)
MCF_WT_resSig_padj <- MCF_WT_res[which(MCF_WT_res$padj < 0.05 ), ]
#rowcounts 22603 (nrows)
MCF10A_H1047R_resSig_padj <- MCF10A_H1047R_res[which(MCF10A_H1047R_res$padj < 0.05 ), ]
#rowcounts 17635 (nrows)
MCF10A_Mut_resSig_padj <- MCF10A_Mut_res[which(MCF10A_Mut_res$padj < 0.05 ), ]
#rowcounts 16924 (nrows)

MCF7_resSig_l2fc <- MCF7_res[which(MCF7_res$log2FoldChange < 0.05 ), ]
#rowcounts 22526 (nrows)
MCF10A_E545K_resSig_l2fc <- MCF10A_E545K_res[which(MCF10A_E545K_res$log2FoldChange < 0.05 ), ]
#rowcounts 24292 (nrows)
MCF_E545K_resSig_l2fc <- MCF_E545K_res[which(MCF_E545K_res$log2FoldChange < 0.05 ), ]
#rowcounts 22788 (nrows)
MCF_WT_resSig_l2fc <- MCF_WT_res[which(MCF_WT_res$log2FoldChange < 0.05 ), ]
#rowcounts 23963 (nrows)
MCF10A_H1047R_resSig_l2fc <- MCF10A_H1047R_res[which(MCF10A_H1047R_res$log2FoldChange < 0.05 ), ]
#rowcounts 25506 (nrows)
MCF10A_Mut_resSig_l2fc <- MCF10A_Mut_res[which(MCF10A_Mut_res$log2FoldChange < 0.05 ), ]
#rowcounts 25868 (nrows)

#using normalization factors (size factors) to account for differences in RNA library depth
MCF_total_dds <- estimateSizeFactors(MCF_total_dds)

#generating a transformed version of count data from samples using vsd or rld, which is useful for downstream analyses -- e.g., for visualization or clustering (PCA, pheatmap)
#varianceStabilizingTransformation = VST
#vst tansformation is less sensitive to outliers, best for medium to large datasets
#This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size
#setting blind to FALSE means the transformation will NOT be blind to the sample information provided by the design formula
MCF_total_vst <- vst(MCF_total_dds, blind = FALSE)
head(assay(MCF_total_vst), 3)
#the assay function extracts matrix of normalized values

#works well on small datasets and outperforms vst when there is a wide range of sequencing depth (not necessarily us but the sample number is --> will use this)
#regularizedlogarithm = rlog
#rlog also produces transformed data on the log2 scale like the vst function
#rlog-transformed data are approximately homoskedastic, so the variance does not depenend on the mean
MCF_total_rlog <- rlog(MCF_total_dds, blind = FALSE)
head(assay(MCF_total_rlog), 3)

#calculate the Euclidean distance between samples, use VST to confirm ~equal contribution from all genes
#Euclidean distance is the distance between two points in a two-dimensional space
#Heatmap shows the Euclidean distance between samples as calculated from the VST of the count data
#A heatmap of this distance matrix gives us an overview of similarities and dissimilarities between samples
#pheatmap is a function of the pheatmap package
dev.off()
sampleDists <- dist(t(assay(MCF_total_vst)))
#extracts the rlog matrix from the object

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- NULL
colnames(sampleDistMatrix) <- paste(MCF_total_samples)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists)

#creating a PCA plot of all replicates = how well do our replicates cluster together? does our experimental condition represent the major source of variation in our data?
#PCA plots are used to evalaute the relations among samples
#In other words, PCA summarizes variance in our dataset, or how our samples vary from each other
#PC1 vs PC2: PC1 is x-axis and describes the direction which separates the data points the most (i.e., largest spread in the data); PC2 represents the second most amount of varation in the data
#The amount of the total variance which is contained in the direction is printed in the axis label
#plotPCA is a function of the DESeq2 package
#88% of the differences among our samples can be grouped together in a way that defines that axis of separation and that other 10% compromises another 'grouping' of separations that fall into the other category
#using vst transformation of the normalized counts as it moderates the variance across the mean, improving the clustering
#by default the plotPCA() function only returns the values for PC1 and PC2
plotPCA(MCF_total_rlog, intgroup = "MCF_total_treatment",)

#plotting histograms of log2fold change and padj
#customize the PCA plot using the ggplot function
MCF7_resSig_padj_df <- as.data.frame(MCF7_resSig_padj)
ggplot(MCF7_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
MCF10A_E545K_resSig_padj_df <- as.data.frame(MCF10A_E545K_resSig_padj)
ggplot(MCF10A_E545K_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
MCF_E545K_resSig_padj_df <- as.data.frame(MCF_E545K_resSig_padj)
ggplot(MCF_E545K_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
MCF_WT_resSig_padj_df <- as.data.frame(MCF_WT_resSig_padj)
ggplot(MCF_WT_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
MCF10A_H1047R_resSig_padj_df <- as.data.frame(MCF10A_H1047R_resSig_padj)
ggplot(MCF10A_H1047R_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
MCF10A_Mut_resSig_padj_df <- as.data.frame(MCF10A_Mut_resSig_padj)
ggplot(MCF10A_Mut_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")

#annotate and export results
#Entrez Gene IDs ("eg") is primary key
library("AnnotationDbi")
library("org.Hs.eg.db")
#add gene symbols and Entrez ID to the res file since it only contains Ensemble gene IDs
#mapIds function adds individual columns to results table
#call mapIds twice in order to add gene symbol and Entez ID
ens.str_7 <- substr(rownames(MCF7_resSig_padj), 1, 15)
ens.str_10_E545K <- substr(rownames(MCF10A_E545K_resSig_padj), 1, 15)
ens.str_E545K <- substr(rownames(MCF_E545K_resSig_padj), 1, 15)
ens.str_WT <- substr(rownames(MCF_WT_resSig_padj), 1, 15)
ens.str_10_H1047R <- substr(rownames(MCF10A_H1047R_resSig_padj), 1, 15)
ens.str_10_Mut <- substr(rownames(MCF10A_Mut_resSig_padj), 1, 15)

#Add Gene Symbols
MCF7_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                   keys=ens.str_7,
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")
MCF10A_E545K_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                    keys=ens.str_10_E545K,
                                    column="SYMBOL",
                                    keytype="ENSEMBL",
                                    multiVals="first")
MCF_E545K_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                        keys=ens.str_E545K,
                                        column="SYMBOL",
                                        keytype="ENSEMBL",
                                        multiVals="first")
MCF_WT_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                    keys=ens.str_WT,
                                    column="SYMBOL",
                                    keytype="ENSEMBL",
                                    multiVals="first")

MCF10A_H1047R_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                    keys=ens.str_10_H1047R,
                                    column="SYMBOL",
                                    keytype="ENSEMBL",
                                    multiVals="first")

MCF10A_Mut_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                    keys=ens.str_10_Mut,
                                    column="SYMBOL",
                                    keytype="ENSEMBL",
                                    multiVals="first")

#Add Entrez IDs
MCF7_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                                   keys=ens.str_7,
                                   column="ENTREZID",
                                   keytype="ENSEMBL",
                                   multiVals="first")
MCF10A_E545K_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                                    keys=ens.str_10_E545K,
                                    column="ENTREZID",
                                    keytype="ENSEMBL",
                                    multiVals="first")
MCF_E545K_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                                        keys=ens.str_E545K,
                                        column="ENTREZID",
                                        keytype="ENSEMBL",
                                        multiVals="first")

MCF_WT_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                               keys=ens.str_WT,
                               column="ENTREZID",
                               keytype="ENSEMBL",
                               multiVals="first")

MCF10A_H1047R_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                               keys=ens.str_10_H1047R,
                               column="ENTREZID",
                               keytype="ENSEMBL",
                               multiVals="first")

MCF10A_Mut_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                               keys=ens.str_10_Mut,
                               column="ENTREZID",
                               keytype="ENSEMBL",
                               multiVals="first")

#ordering the new resSig files according to their pvalue
MCF7_resSig_padj_ordered <- MCF7_resSig_padj[order(MCF7_resSig_padj$pvalue),]
MCF7_resSig_padj_ordered_df <- as.data.frame(MCF7_resSig_padj_ordered)

MCF10A_E545K_resSig_padj_ordered <- MCF10A_E545K_resSig_padj[order(MCF10A_E545K_resSig_padj$pvalue),]
MCF10A_E545K_resSig_padj_ordered_df <- as.data.frame(MCF10A_E545K_resSig_padj_ordered)

MCF_E545K_resSig_padj_ordered <- MCF_E545K_resSig_padj[order(MCF_E545K_resSig_padj$pvalue),]
MCF_E545K_resSig_padj_ordered_df <- as.data.frame(MCF_E545K_resSig_padj_ordered)

MCF_WT_resSig_padj_ordered <- MCF_WT_resSig_padj[order(MCF_WT_resSig_padj$pvalue),]
MCF_WT_resSig_padj_ordered_df <- as.data.frame(MCF_WT_resSig_padj_ordered)

MCF10A_H1047R_resSig_padj_ordered <- MCF10A_H1047R_resSig_padj[order(MCF10A_H1047R_resSig_padj$pvalue),]
MCF10A_H1047R_resSig_padj_ordered_df <- as.data.frame(MCF10A_H1047R_resSig_padj_ordered)

MCF10A_Mut_resSig_padj_ordered <- MCF10A_Mut_resSig_padj[order(MCF10A_Mut_resSig_padj$pvalue),]
MCF10A_Mut_resSig_padj_ordered_df <- as.data.frame(MCF10A_Mut_resSig_padj_ordered)

#export these results in a CSV file, which you can load with a spreadsheet program such as Excel
write.csv(MCF7_resSig_padj_ordered_df, file = "MCF7_correctedvparental.csv")
write.csv(MCF10A_E545K_resSig_padj_ordered_df, file = "MCF10_E545Kvparental.csv")
write.csv(MCF_E545K_resSig_padj_ordered_df, file = "MCF10AE545KvMCF7parental.csv")
write.csv(MCF_WT_resSig_padj_ordered_df, file = "MCF10AparentalvMCF7corrected.csv")
write.csv(MCF10A_H1047R_resSig_padj_ordered_df, file = "MCF10A_H1047Rvparental.csv")
write.csv(MCF10A_Mut_resSig_padj_ordered_df, file = "MCF10A_H1047RvE545K.csv")




###Making a MCF10A only comparison table to create a more relevant PCA/heatmap#################################################################
###############################################################################################################################################
###############################################################################################################################################
MCF10A_counts <- total_counts[c("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")]
MCF_10A_samples <- c("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
MCF_10A_treatment<- c("10AE545K", "10AE545K", "10AE545K", "10AParental", "10AParental", "10AParental", "10AH1047R", "10AH1047R", "10AH1047R")
MCF_10A_treatment <- as.factor(MCF_10A_treatment)
MCF_10A_coldata <- data.frame(MCF_10A_treatment)
row.names(MCF_10A_coldata) <- MCF_10A_samples
all(rownames(MCF_10A_coldata) == colnames(MCF10A_counts))
MCF_10A_dds <- DESeqDataSetFromMatrix(countData= MCF10A_counts,
                                        colData= MCF_10A_coldata,
                                        design= ~ MCF_10A_treatment)
MCF_10A_counts_filtered <- rowSums(counts(MCF_10A_dds)) > 10
MCF_10A_dds <- MCF_10A_dds[MCF_10A_counts_filtered,]
MCF_10A_dds <- estimateSizeFactors(MCF_10A_dds)
MCF_10A_vst <- vst(MCF_10A_dds, blind = FALSE)
head(assay(MCF_10A_vst), 3)
MCF_10A_rlog <- rlog(MCF_10A_dds, blind = FALSE)
head(assay(MCF_10A_rlog), 3)
sampleDists_10 <- dist(t(assay(MCF_10A_vst)))
sampleDistMatrix_10 <- as.matrix( sampleDists_10 )
rownames(sampleDistMatrix_10) <- NULL
colnames(sampleDistMatrix_10) <- paste(MCF_10A_samples)
pheatmap(sampleDistMatrix_10,
         clustering_distance_rows = sampleDists_10,
         clustering_distance_cols = sampleDists_10)
plotPCA(MCF_10A_rlog, intgroup = "MCF_10A_treatment",)
#################################################################
#################################################################
#################################################################