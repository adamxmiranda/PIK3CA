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

#read in the feature counts file
#read.delim is used to read in delimited text files, where data is organized in a data matrix with rows representing cases and columns representing variables
#select() allows you to zoom in on columns that are actually of interest to you in the feature counts file using dyplyr
MCF_dual_counts <- read.delim(file = "~/Desktop/RNAseq_dual/featureCounts_all_samples_MCF.txt", skip = 1,
                           col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "MCF10A_mutant_rep1", "MCF10A_mutant_rep2", "MCF10A_mutant_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3"),
                           row.names = "geneID")  %>%
  dplyr::select("MCF10A_mutant_rep1", "MCF10A_mutant_rep2", "MCF10A_mutant_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3")

#creating a coldata table for the DESeq matrix, which is just a table with metadata on the count table's columns; colData = column metadata table
#colData = metadata = the "columns" of the count table
#colData is a data.frame that can contain all the variables you know about your samples, such as the experimental condition, the type and date of sequencing and so on
MCF_dual_samples <- c("MCF10A_mutant_rep1", "MCF10A_mutant_rep2", "MCF10A_mutant_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3")
MCF_dual_treatment<- c("10AMutant", "10AMutant", "10AMutant", "10AParental", "10AParental", "10AParental", "7Corrected", "7Corrected", "7Corrected", "7Parental","7Parental", "7Parental")

#values need to be in factor format for DESeq
MCF_dual_treatment <- as.factor(MCF_dual_treatment)

#character variables passed to data.frame are converted to factor columns (i.e., CRISPRi treatment)
#row names are CRISPRi samples from colData data frame
MCF_dual_coldata <- data.frame(MCF_dual_treatment)
row.names(MCF_dual_coldata) <- MCF_dual_samples
all(rownames(MCF_dual_coldata) == colnames(MCF_dual_counts))

#creating DESeq data object from the matrix of counts and the metadata table (colData) using DESeqDataSetFromMatrix input function
#testing for the effect of treatment (e.g., Control CRISPRi versus CRISPRi targeting ZNF384) with design variable CRISPRi_treatment
#countData = siRNA_counts = a table with the full read counts
#to avoid model matrix is not a full rank error can only consider the treatment variable in the design
#design formula accounts for source of variation (i.e., CRISPRi_treatment)
MCF_dual_dds <- DESeqDataSetFromMatrix(countData= MCF_dual_counts,
                                   colData= MCF_dual_coldata,
                                   design= ~ MCF_dual_treatment)

#To run the differential expression analysis, we use a single call function DESeq()
#This function prints out a message for various steps: estimating pre-size factors, dispersons etc...
#has every comparison that can be drawn
MCF_dual_dds <- DESeq(MCF_dual_dds)
View(counts(MCF_dual_dds))

MCF_dual_dds_df <- counts(MCF_dual_dds)
write.csv(MCF_dual_dds_df, file = "MCF_expression_BRAD.csv")
#to see what comparisons were made; resultsNames returns the names of the estimated effects (coefficients) of the model
resultsNames(MCF_dual_dds)

#To normalize the count data, DESeq2 calculates size factors (i.e., normalization factors) for each samples using the median of ratios method, which is automatically done when performing the differential expression analysis -- DESeq()
#Normalization takes into account: sequencing depth, gene length, and RNA composition
#Median of ratios method makes the assumptions that not ALL genes are differentially expressed
#Check the size factors; each sample has a unique normalization (size) factor
sizeFactors(MCF_dual_dds)

#Total number of raw counts per sample
#How do the numbers correlate with size factor?
colSums(counts(MCF_dual_dds))

#Total number of normalized counts per sample
#Note: DESeq2 doesn't actually used normalized counts; need raw counts as input (i.e., featureCounts txt file)
colSums(counts(MCF_dual_dds, normalized=T))
normalized_MCF_dual_counts <- counts(MCF_dual_dds, normalized=T)

#pre-filtering the counts so to eliminates genes that had no reads or less than 10 reads map to them
#kick out genes that have a total read count less than 10 across all samples, which leaves 34982 genes out of the total count
#pre-filtering low gene counts reduces the memory size of the dds data object (CRISPRi_dds) and increases the speed of the transformation and testing functions within DESeq2
MCF_dual_counts_filtered <- rowSums(counts(MCF_dual_dds)) > 10
MCF_dual_dds <- MCF_dual_dds[MCF_dual_counts_filtered,]
#rowcount 41879

#save dds file
save(MCF_dual_dds, file = "MCF_dual_dds.Rdata")

#independent filtering is performed automatically with the results function
#outputting the filtered results, which extracts a results table with log2 fold changes, p values and adjusted p values
#in a previous step, we filtered out the genes that that zero counts in all samples or less than 10 reads across
#now with results() we  remove genes with an extreme count outlier and/or a low mean normalized count
#reference level = Control_CRISPRi
#e.g., log2(Parental/Corrected CRISPRi)
MCF_7_res <- results(MCF_dual_dds, contrast = c("MCF_dual_treatment", "7Corrected", "7Parental"))
MCF_10_res <- results(MCF_dual_dds, contrast = c("MCF_dual_treatment", "10AMutant", "10AParental"))
MCF_Mutant_res <- results(MCF_dual_dds, contrast = c("MCF_dual_treatment", "10AMutant", "7Parental"))
MCF_WT_res <- results(MCF_dual_dds, contrast = c("MCF_dual_treatment", "10AParental", "7Corrected"))
#variable name, numerator, denominator

MCF_7_res_df <- as.data.frame(MCF_7_res)
MCF_10_res_df <- as.data.frame(MCF_10_res)
MCF_Mutant_res_df <- as.data.frame(MCF_Mutant_res)
MCF_WT_res_df <- as.data.frame(MCF_WT_res)
#data.frame allows you to use CRISPRi results for other R packages

diff_7 <- rownames_to_column(MCF_7_res_df)
diff_10 <- rownames_to_column(MCF_10_res_df)
diff_Mutant <- rownames_to_column(MCF_Mutant_res_df)
diff_WT <- rownames_to_column(MCF_WT_res_df)
#row names become a column

mcols(MCF_7_res, use.names=TRUE)
mcols(MCF_10_res, use.names=TRUE)
mcols(MCF_Mutant_res, use.names=TRUE)
mcols(MCF_WT_res, use.names=TRUE)
#carries metadata with information on the meaning of the columns
#dif_2 <- diff %>% dplyr::filter(str_detect(diff$rowname, pattern = "ENSG00000126746")) --> ENSEMBLE code for ZNF384

#filtering the results by their padj value
#if we consider a fraction of 5% false positives acceptable among differentially expressed genes, we can consider all genes with an adjusted p value below 5%=0.05 as significant; how many such genes are there?
#A separate adjusted p value is computed for each comparsion in a family of comparsions and each comparision has a unique adjusted P value
#P adj is the smallest family wise significance level at which a particular comparsion will be declared statistically significant as part of multiple comparison testing
MCF_7_resSig_padj <- MCF_7_res[which(MCF_7_res$padj < 0.05 ), ]
#rowcounts 14365 (nrows)
MCF_10_resSig_padj <- MCF_10_res[which(MCF_10_res$padj < 0.05 ), ]
#rowcounts 12209 (nrows)
MCF_Mutant_resSig_padj <- MCF_Mutant_res[which(MCF_Mutant_res$padj < 0.05 ), ]
#rowcounts 22710 (nrows)
MCF_WT_resSig_padj <- MCF_WT_res[which(MCF_WT_res$padj < 0.05 ), ]
#rowcounts 23292 (nrows)

MCF_7_resSig_l2fc <- MCF_7_res[which(MCF_7_res$log2FoldChange >=0), ]
#rowcounts 22338 (nrows)
MCF_10_resSig_l2fc <- MCF_10_res[which(MCF_10_res$log2FoldChange >=0), ]
#rowcounts 21083 (nrows)
MCF_Mutant_resSig_l2fc <- MCF_Mutant_res[which(MCF_Mutant_res$log2FoldChange >=0), ]
#rowcounts 21119 (nrows)
MCF_WT_resSig_l2fc <- MCF_WT_res[which(MCF_WT_res$log2FoldChange >=0), ]
#rowcounts 20009 (nrows)

#using normalization factors (size factors) to account for differences in RNA library depth
MCF_dual_dds <- estimateSizeFactors(MCF_dual_dds)

#generating a transformed version of count data from samples using vsd or rld, which is useful for downstream analyses -- e.g., for visualization or clustering (PCA, pheatmap)
#varianceStabilizingTransformation = VST
#vst tansformation is less sensitive to outliers, best for medium to large datasets
#This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size
#setting blind to FALSE means the transformation will NOT be blind to the sample information provided by the design formula
MCF_dual_vst <- vst(MCF_dual_dds, blind = FALSE)
head(assay(MCF_dual_vst), 3)
#the assay function extracts matrix of normalized values

#works well on small datasets and outperforms vst when there is a wide range of sequencing depth (not necessarily us but the sample number is --> will use this)
#regularizedlogarithm = rlog
#rlog also produces transformed data on the log2 scale like the vst function
#rlog-transformed data are approximately homoskedastic, so the variance does not depenend on the mean
MCF_dual_rlog <- rlog(MCF_dual_dds, blind = FALSE)
head(assay(MCF_dual_rlog), 3)

#calculate the Euclidean distance between samples, use VST to confirm ~equal contribution from all genes
#Euclidean distance is the distance between two points in a two-dimensional space
#Heatmap shows the Euclidean distance between samples as calculated from the VST of the count data
#A heatmap of this distance matrix gives us an overview of similarities and dissimilarities between samples
#pheatmap is a function of the pheatmap package
dev.off()
sampleDists <- dist(t(assay(MCF_dual_vst)))
#extracts the rlog matrix from the object

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- NULL
colnames(sampleDistMatrix) <- paste(MCF_dual_samples)
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
plotPCA(MCF_dual_rlog, intgroup = "MCF_dual_treatment",)

#plotting histograms of log2fold change and padj
#customize the PCA plot using the ggplot function
MCF_7_resSig_padj_df <- as.data.frame(MCF_7_resSig_padj)
ggplot(MCF_7_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
MCF_10_resSig_padj_df <- as.data.frame(MCF_10_resSig_padj)
ggplot(MCF_10_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
MCF_Mutant_resSig_padj_df <- as.data.frame(MCF_Mutant_resSig_padj)
ggplot(MCF_Mutant_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
MCF_WT_resSig_padj_df <- as.data.frame(MCF_WT_resSig_padj)
ggplot(MCF_WT_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")

#lets try clustering shall we
MCF_dual_rlog_mat <- as.matrix(assay(MCF_dual_rlog))
gap_stat <- clusGap(MCF_dual_rlog_mat, FUN = kmeans, nstart = 30, K.max = 25, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")


###Make Venn Diagrams!!!
library(VennDiagram)
library(RColorBrewer)
MCF_7_genes <- rownames(MCF_7_resSig_padj_df)
MCF_10_genes <- rownames(MCF_10_resSig_padj_df)
MCF_Mutant_genes <- rownames(MCF_Mutant_resSig_padj_df)
MCF_WT_genes <- rownames(MCF_WT_resSig_padj_df)
colors <- brewer.pal(4, "Paired")
venn.diagram(
  x= list(MCF_7_genes, MCF_10_genes, MCF_Mutant_genes, MCF_WT_genes),
  category.names = c("MCF7 Parental v Corrected" , "MCF10A Parental v Mutant" , "MCF7 Parental v MCF10A Mutant", "MCF10A Parental v MCF7 Corrected"),
  filename = "Venn.png",
  height = 5000,
  width = 5000,
  output= FALSE,
  col = colors,
  fill= alpha(colors, 0.5),
  cat.fontface= "bold")

#annotate and export results
#Entrez Gene IDs ("eg") is primary key
library("AnnotationDbi")
library("org.Hs.eg.db")
#add gene symbols and Entrez ID to the res file since it only contains Ensemble gene IDs
#mapIds function adds individual columns to results table
#call mapIds twice in order to add gene symbol and Entez ID
ens.str_7 <- substr(rownames(MCF_7_resSig_padj), 1, 15)
ens.str_10 <- substr(rownames(MCF_10_resSig_padj), 1, 15)
ens.str_Mut <- substr(rownames(MCF_Mutant_resSig_padj), 1, 15)
ens.str_WT <- substr(rownames(MCF_WT_resSig_padj), 1, 15)


MCF_7_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                               keys=ens.str_7,
                               column="SYMBOL",
                               keytype="ENSEMBL",
                               multiVals="first")
MCF_10_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                   keys=ens.str_10,
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")
MCF_Mutant_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                   keys=ens.str_Mut,
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")
MCF_WT_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                   keys=ens.str_WT,
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")

MCF_7_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                               keys=ens.str_7,
                               column="ENTREZID",
                               keytype="ENSEMBL",
                               multiVals="first")
MCF_10_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                              keys=ens.str_10,
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")
MCF_Mutant_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                              keys=ens.str_Mut,
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")
MCF_WT_resSig$entrez <- mapIds(org.Hs.eg.db,
                              keys=ens.str_WT,
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")


#ordering the new resSig files according to their pvalue
MCF_7_resSig_ordered <- MCF_7_resSig_padj[order(MCF_7_resSig_padj$pvalue),]
MCF_7_resSig_ordered_df <- as.data.frame(MCF_7_resSig_ordered)

MCF_10_resSig_ordered <- MCF_10_resSig_padj[order(MCF_10_resSig_padj$pvalue),]
MCF_10_resSig_ordered_df <- as.data.frame(MCF_10_resSig_ordered)

MCF_Mutant_resSig_ordered <- MCF_Mutant_resSig_padj[order(MCF_Mutant_resSig_padj$pvalue),]
MCF_Mutant_resSig_ordered_df <- as.data.frame(MCF_Mutant_resSig_ordered)

MCF_WT_resSig_ordered <- MCF_WT_resSig_padj[order(MCF_WT_resSig_padj$pvalue),]
MCF_WT_resSig_ordered_df <- as.data.frame(MCF_WT_resSig_ordered)
#export these results in a CSV file, which you can load with a spreadsheet program such as Excel
write.csv(MCF_7_resSig_ordered_df, file = "MCF7_correctedvparental.csv")
write.csv(MCF_10_resSig_ordered_df, file = "MCF10_mutantvparental.csv")
write.csv(MCF_Mutant_resSig_ordered_df, file = "MCF10AmutantvMCF7parental.csv")
write.csv(MCF_WT_resSig_ordered_df, file = "MCF10AparentalvMCF7corrected.csv")
