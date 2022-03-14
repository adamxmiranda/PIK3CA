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
library("RColorBrewer")
setwd('~/Desktop/RNAseq_Dual')
#read in the feature counts file
#read.delim is used to read in delimited text files, where data is organized in a data matrix with rows representing cases and columns representing variables
#select() allows you to zoom in on columns that are actually of interest to you in the feature counts file using dyplyr
MCF_dual_E545K_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_MCF.txt", skip = 1,
                              col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3"),
                              row.names = "geneID")  %>%
  dplyr::select("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3", "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3", "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3", "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3")

MCF10A_H1047R_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_MCF10A_H1047R.txt", skip = 1,
                              col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3"),
                              row.names = "geneID")  %>%
  dplyr::select("MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
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
T47D_counts <- read.delim(file = "~/Desktop/Total_RNA/featureCounts_all_samples_T47D.txt", skip = 1,
                                 col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "T47D_rep1", "T47D_rep2", "T47D_rep3"),
                                 row.names = "geneID") %>%
  dplyr::select("T47D_rep1", "T47D_rep2", "T47D_rep3")

#combine dataframes with cbind
total_counts <- cbind(MCF_dual_E545K_counts, MCF10A_H1047R_counts, HTert_WT_counts, HTert_ex9_counts, HTert_ex20_counts, T47D_counts)
#creating a coldata table for the DESeq matrix, which is just a table with metadata on the count table's columns; colData = column metadata table
#colData = metadata = the "columns" of the count table
#colData is a data.frame that can contain all the variables you know about your samples, such as the experimental condition, the type and date of sequencing and so on
total_samples <- c("MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                       "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                       "MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                       "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                       "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3",
                       "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                       "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                       "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
                       "T47D_rep1", "T47D_rep2", "T47D_rep3")
total_treatment <- c("ex9", "ex9", "ex9",
                        "WT", "WT", "WT",
                        "WT", "WT", "WT",
                        "ex9", "ex9", "ex9",
                        "ex20", "ex20", "ex20",
                        "WT", "WT", "WT",
                        "ex9", "ex9", "ex9",
                        "ex20", "ex20", "ex20",
                        "ex20", "ex20", "ex20")

total_treatment_v2 <- c("ex9_MCF10A", "ex9_MCF10A", "ex9_MCF10A",
                    "WT_MCF10A", "WT_MCF10A", "WT_MCF10A",
                    "WT_MCF7", "WT_MCF7", "WT_MCF7",
                    "ex9_MCF7", "ex9_MCF7", "ex9_MCF7",
                    "ex20_MCF10A", "ex20_MCF10A", "ex20_MCF10A",
                    "WT_HTert", "WT_HTert", "WT_HTert",
                    "ex9_HTert", "ex9_HTert", "ex9_HTert",
                    "ex20_HTert", "ex20_HTert", "ex20_HTert",
                    "ex20_T47D", "ex20_T47D", "ex20_T47D")


#values need to be in factor format for DESeq
total_treatment <- as.factor(total_treatment)
total_treatment_v2 <- as.factor(total_treatment_v2)
total_treatment_v3 <- rbind(total_treatment, total_treatment_v2)
#t turns rows into columns 
total_treatment_v3 <- t(total_treatment_v3)

#character variables passed to data.frame are converted to factor columns (i.e., CRISPRi treatment)
#row names are CRISPRi samples from colData data frame
total_coldata <- data.frame(total_treatment)
total_coldata_v2 <- data.frame(total_treatment_v2)
total_coldata_v3 <- data.frame(total_treatment_v3)
row.names(total_coldata) <- total_samples
row.names(total_coldata_v2) <- total_samples
row.names(total_coldata_v3) <- total_samples
all(rownames(total_coldata)) == colnames(total_counts)) #TRUE
all(rownames(total_coldata_v2) == colnames(total_counts)) #TRUE 
all(rownames(total_coldata_v3) == colnames(total_counts)) #TRUE 

#creating DESeq data object from the matrix of counts and the metadata table (colData) using DESeqDataSetFromMatrix input function
#testing for the effect of treatment (e.g., Control CRISPRi versus CRISPRi targeting ZNF384) with design variable CRISPRi_treatment
#countData = siRNA_counts = a table with the full read counts
#to avoid model matrix is not a full rank error can only consider the treatment variable in the design
#design formula accounts for source of variation (i.e., CRISPRi_treatment)
total_dds <- DESeqDataSetFromMatrix(countData= total_counts,
                                       colData= total_coldata,
                                       design= ~ total_treatment)
total_dds_v2 <- DESeqDataSetFromMatrix(countData= total_counts,
                                    colData= total_coldata_v3,
                                    design= ~ total_treatment) 


#collapse by genotype: WT, ex9, and ex20 
#collapseReplicates assists in combining the counts from technical replicates into single columns of the count matrix
total_dds_collapsed <-collapseReplicates(total_dds,
                                         groupby = total_dds$total_treatment) 
colData(total_dds_collapsed)
colnames(total_dds_collapsed)

#collapse by genotype within each cell lineage: ex20_HTert ex20_MCF10A ex20_T47D ex9_HTert ex9_MCF10A ex9_MCF7 WT_HTert WT_MCF10A WT_MCF7
total_dds_collapsed_v2 <- collapseReplicates(total_dds_v2,
                                         groupby = total_dds_v2$total_treatment_v2,
                                         run = total_coldata_v2$total_treatment)
View( as.data.frame( colData(total_dds_collapsed_v2)[ ,c("total_treatment","runsCollapsed") ] ))

as.data.frame( colData(total_dds_collapsed_v2) )

###Additional Pre-filtering step for low count genes
filt <- rowSums( counts(total_dds, normalized=FALSE) >= 20 ) >= 14
filt_v2 <- rowSums( counts(total_dds_collapsed_v2, normalized=FALSE) >= 20 ) >=5
##use 10
filt_v3 <- rowSums( counts(total_dds_collapsed_v2, normalized=FALSE) >= 10 ) >=5


total_dds_f <- total_dds[filt,]
total_dds_f 
total_dds_f <- DESeq(total_dds_f)

total_dds_collapsed_f <- total_dds_collapsed_v2[filt_v2,]
total_dds_collapsed_f_v3 <- total_dds_collapsed_v2[filt_v3,]
nrow(total_dds_collapsed_f) #26235
nrow(total_dds_collapsed_f_v3) #30262


#can't collapse by genotype since in order to do statistical analysis, you need replicates to estimate the biological variability
total_dds_collapsed_f <- DESeq(total_dds_collapsed_f)
total_dds_collapsed_f_v3_DESeq <- DESeq(total_dds_collapsed_f_v3)
collapsed_res <- results(total_dds_collapsed_f)


normalized_total_counts_f <- counts(total_dds_f, normalized=T)
View(normalized_total_counts_f)

normalized_total_counts_f_nonames <- cbind(rownames(normalized_total_counts_f), data.frame(normalized_total_counts_f, row.names=NULL))
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
total_dds <- DESeq(total_dds_f)
total_dds_collapsed_DESeq <- DESeq(total_dds_collapsed_f)

#rows with 0s were removed 
x <- as.data.frame(counts(total_dds_collapsed_f_v3_DESeq))
nrow(x)

#to see what comparisons were made; resultsNames returns the names of the estimated effects (coefficients) of the model
resultsNames(total_dds)
resultsNames(total_dds_collapsed_f_v3_DESeq)
#To normalize the count data, DESeq2 calculates size factors (i.e., normalization factors) for each samples using the median of ratios method, which is automatically done when performing the differential expression analysis -- DESeq()
#Normalization takes into account: sequencing depth, gene length, and RNA composition
#Median of ratios method makes the assumptions that not ALL genes are differentially expressed
#Check the size factors; each sample has a unique normalization (size) factor
sizeFactors(total_dds)
sizeFactors(total_dds_collapsed_DESeq)
#Total number of raw counts per sample
#How do the numbers correlate with size factor?
colSums(counts(total_dds))
#Total number of normalized counts per sample
#Note: DESeq2 doesn't actually used normalized counts; need raw counts as input (i.e., featureCounts txt file)
colSums(counts(total_dds, normalized=T))
normalized_total_counts <- counts(total_dds, normalized=T)

#save dds file
save(total_dds, file = "total_dds.Rdata")
save(total_dds_collapsed_f_v3_DESeq, file = "total_dds_collapsed_f_v3_DESeq.Rdata")

#independent filtering is performed automatically with the results function
#outputting the filtered results, which extracts a results table with log2 fold changes, p values and adjusted p values
#in a previous step, we filtered out the genes that have zero counts in all samples or less than 10 reads across
#now with results() we  remove genes with an extreme count outlier and/or a low mean normalized count
#reference level = Control_CRISPRi
#e.g., log2(Parental/Corrected CRISPRi)
#this is where the comparisons are made, list all of the comparisons that are relevant to the data
#By default the argument alpha is set to 0.1.
WTvex9_res <- results(total_dds_collapsed_f_v3_DESeq, contrast = c("total_treatment", "WT", "ex9"))
WTvex20_res <- results(total_dds_collapsed_f_v3_DESeq, contrast = c("total_treatment", "WT", "ex20"))
ex9vex20_res <- results(total_dds_collapsed_f_v3_DESeq, contrast = c("total_treatment", "ex9", "ex20"))

x <- plotMA(WTvex9_res, ylim=c(-2,2))
y <- plotMA(WTvex20_res, ylim=c(-2,2))
z <- plotMA(ex9vex20_res, ylim=c(-2,2))

list <- c(WTvex9_res, WTvex20_res, ex9vex20_res)

resOrdered <- WTvex9_res[order(WTvex9_res$padj),]

padj_filter <- function(inDF){
  print(sum(inDF$padj < 1.0, na.rm=TRUE))
}
padj_filter(WTvex9_res)

#variable name, numerator, denominator
WTvex9_res_df <- as.data.frame(WTvex9_res)
WTvex20_res_df <- as.data.frame(WTvex20_res)
ex9vex20_res_df <- as.data.frame(ex9vex20_res)

#data.frame allows you to use CRISPRi results for other R packages
diff_WTvex9 <- rownames_to_column(WTvex9_res_df)
diff_WTvex20 <- rownames_to_column(WTvex20_res_df)
diff_ex9vex20 <- rownames_to_column(ex9vex20_res_df)

#row names become a column
mcols(WTvex9_res, use.names=TRUE)
mcols(WTvex20_res, use.names=TRUE)
mcols(ex9vex20_res, use.names=TRUE)

#carries metadata with information on the meaning of the columns
#dif_2 <- diff %>% dplyr::filter(str_detect(diff$rowname, pattern = "ENSG00000126746")) --> ENSEMBLE code for ZNF384

WTvex9_resSig_padj <- WTvex9_res[which(WTvex9_res$pvalue < 0.1 ), ]
#rowcounts 447 (nrows)

WTvex20_resSig_padj <- WTvex20_res[which(WTvex20_res$padj < 0.1 ), ]
#rowcounts 1944 (nrows)

ex9vex20_resSig_padj <- ex9vex20_res[which(ex9vex20_res$padj < 0.1 ), ]
#rowcounts 2388 (nrows)

#using normalization factors (size factors) to account for differences in RNA library depth
total_dds <- estimateSizeFactors(total_dds)

#generating a transformed version of count data from samples using vsd or rld, which is useful for downstream analyses -- e.g., for visualization or clustering (PCA, pheatmap)
#varianceStabilizingTransformation = VST
#vst tansformation is less sensitive to outliers, best for medium to large datasets
#This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size
#setting blind to FALSE means the transformation will NOT be blind to the sample information provided by the design formula
total_vst <- vst(total_dds, blind = FALSE)
head(assay(total_vst), 3)
#the assay function extracts matrix of normalized values

#works well on small datasets and outperforms vst when there is a wide range of sequencing depth (not necessarily us but the sample number is --> will use this)
#regularizedlogarithm = rlog
#rlog also produces transformed data on the log2 scale like the vst function
#rlog-transformed data are approximately homoskedastic, so the variance does not depenend on the mean
total_rlog <- rlog(total_dds, blind = FALSE)
head(assay(total_rlog), 3)

#calculate the Euclidean distance between samples, use VST to confirm ~equal contribution from all genes
#Euclidean distance is the distance between two points in a two-dimensional space
#Heatmap shows the Euclidean distance between samples as calculated from the VST of the count data
#A heatmap of this distance matrix gives us an overview of similarities and dissimilarities between samples
#pheatmap is a function of the pheatmap package
dev.off()
sampleDists <- dist(t(assay(total_vst)))
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
plotPCA(total_rlog, intgroup = "total_treatment",)

#plotting histograms of log2fold change and padj
#customize the PCA plot using the ggplot function
WTvex9_resSig_padj_df <- as.data.frame(WTvex9_resSig_padj)
ggplot(WTvex9_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
WTvex20_resSig_padj_df <- as.data.frame(WTvex20_resSig_padj)
ggplot(WTvex20_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")
ex9vex20_resSig_padj_df <- as.data.frame(ex9vex20_resSig_padj)
ggplot(ex9vex20_resSig_padj_df, aes(x=log2FoldChange)) +
  geom_histogram(binwidth=.1, color="black", fill="gray")

#annotate and export results
#Entrez Gene IDs ("eg") is primary key
library("AnnotationDbi")
library("org.Hs.eg.db")
#add gene symbols and Entrez ID to the res file since it only contains Ensemble gene IDs
#mapIds function adds individual columns to results table
#call mapIds twice in order to add gene symbol and Entez ID
ens.str_WTvex9 <- substr(rownames(WTvex9_resSig_padj), 1, 15)
ens.str_WTvex20 <- substr(rownames(WTvex20_resSig_padj), 1, 15)
ens.str_ex9vex20 <- substr(rownames(ex9vex20_resSig_padj), 1, 15)

#Add Gene Symbols
WTvex9_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                   keys=ens.str_WTvex9,
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")
WTvex20_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                    keys=ens.str_WTvex20,
                                    column="SYMBOL",
                                    keytype="ENSEMBL",
                                    multiVals="first")
ex9vex20_resSig_padj$symbol <- mapIds(org.Hs.eg.db,
                                        keys=ens.str_ex9vex20,
                                        column="SYMBOL",
                                        keytype="ENSEMBL",
                                        multiVals="first")

#Add Entrez IDs
WTvex9_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                                   keys=ens.str_WTvex9,
                                   column="ENTREZID",
                                   keytype="ENSEMBL",
                                   multiVals="first")
WTvex20_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                                    keys=ens.str_WTvex20,
                                    column="ENTREZID",
                                    keytype="ENSEMBL",
                                    multiVals="first")
ex9vex20_resSig_padj$entrez <- mapIds(org.Hs.eg.db,
                                        keys=ens.str_ex9vex20,
                                        column="ENTREZID",
                                        keytype="ENSEMBL",
                                        multiVals="first")


#ordering the new resSig files according to their pvalue
WTvex9_resSig_padj_ordered <- WTvex9_resSig_padj[order(WTvex9_resSig_padj$pvalue),]
WTvex9_resSig_padj_ordered_df <- as.data.frame(WTvex9_resSig_padj_ordered)

WTvex20_resSig_padj_ordered <- WTvex20_resSig_padj[order(WTvex20_resSig_padj$pvalue),]
WTvex20_resSig_padj_ordered_df <- as.data.frame(WTvex20_resSig_padj_ordered)

ex9vex20_resSig_padj_ordered <- ex9vex20_resSig_padj[order(ex9vex20_resSig_padj$pvalue),]
ex9vex20_resSig_padj_ordered_df <- as.data.frame(ex9vex20_resSig_padj_ordered)

#export these results in a CSV file, which you can load with a spreadsheet program such as Excel
write.csv(WTvex9_resSig_padj_ordered_df, file = "WTvex9.csv")
write.csv(WTvex20_resSig_padj_ordered_df, file = "WTvex20.csv")
write.csv(ex9vex20_resSig_padj_ordered_df, file = "ex9vex20.csv")



###Attempt at k-means clustering on significant subset###
#create gene list of differentially expressed genes across all three comparisons
#foldchange and p-value filter
diff_WTvex9_sig <- diff_WTvex9 %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
diff_WTvex9_sig <- diff_WTvex9_sig %>%
  filter(padj < .05)

diff_WTvex20_sig <- diff_WTvex20 %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
diff_WTvex20_sig <- diff_WTvex20_sig %>%
  filter(padj < .05)

diff_ex9vex20_sig <- diff_ex9vex20 %>%
  filter(log2FoldChange < -1.5 | log2FoldChange > 1.5)
diff_ex9vex20_sig <- diff_ex9vex20_sig %>%
  filter(padj < .05)

all_diff_genes <- rbind(diff_WTvex9_sig, diff_WTvex20_sig, diff_ex9vex20_sig)
all_diff_genes <- all_diff_genes[!duplicated(all_diff_genes$rowname),]

#create a subset of normalized reads using this gene list
norm_total_counts_sig_subset <- subset(normalized_total_counts_f_nonames, normalized_total_counts_f_nonames$`rownames(normalized_total_counts_f)` %in% all_diff_genes$rowname)
norm_total_counts_sig_subset <- remove_rownames(norm_total_counts_sig_subset)
norm_total_counts_sig_subset <- as.data.frame(norm_total_counts_sig_subset)
norm_total_counts_sig_subset <- column_to_rownames(norm_total_counts_sig_subset, var = "rownames(normalized_total_counts_f)")

#try a raw counts version
#first create a raw counts matrix with no rownames
total_counts_nonames <- cbind(rownames(total_counts), data.frame(total_counts, row.names=NULL))
raw_total_counts_sig_subset <- subset(total_counts_nonames, total_counts_nonames$`rownames(total_counts)` %in% all_diff_genes$rowname)
raw_total_counts_sig_subset <- remove_rownames(raw_total_counts_sig_subset )
raw_total_counts_sig_subset  <- as.data.frame(raw_total_counts_sig_subset )
raw_total_counts_sig_subset  <- column_to_rownames(raw_total_counts_sig_subset , var = "rownames(total_counts)")

#k-means clustering
#use elbow method to determine the number of clusters present in the data
kmax <- 15
wss <-sapply(1:kmax,
             function(k){kmeans(raw_total_counts_sig_subset, k,
                                nstart = 50, iter.max= 15) $tot.withinss})
wss
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
#looks like 6 clusters
#create heatmap with clustering
p <- pheatmap(norm_total_counts_sig_subset, kmeans_k = 6)
p <- pheatmap(raw_total_counts_sig_subset, kmeans_k = 6)
#try clustering by correlation
p <- pheatmap(norm_total_counts_sig_subset, clustering_distance_rows= "correlation", clustering_distance_cols = "correlation")
p <- pheatmap(raw_total_counts_sig_subset, clustering_distance_rows= "correlation", clustering_distance_cols = "correlation")


#AREG Expression
AREG_exp <- as.data.frame(normalized_total_counts_t$ENSG00000109321.11)
AREG_exp$cell_line <- rownames(normalized_total_counts_t)
AREG_exp$AREG_NormCounts <- AREG_exp$`normalized_total_counts_t$ENSG00000109321.11`
reorder_cells <- c("MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                   "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                   "T47D_rep1", "T47D_rep2", "T47D_rep3",
                   "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                   "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                   "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
                   "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                   "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                   "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
AREG_exp <- AREG_exp %>%
  dplyr :: slice(match(reorder_cells, cell_line))
AREG_exp$cell_line <- factor(AREG_exp$cell_line, levels = AREG_exp$cell_line)
ggplot(AREG_exp, aes(x = cell_line,
                     y = AREG_NormCounts)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() +
  coord_flip()

total_counts_t <- t(total_counts)
total_counts_t <- as.data.frame(total_counts_t)
AREG_exp_raw <- as.data.frame(total_counts_t$ENSG00000109321.11)
AREG_exp_raw$cell_line <- rownames(total_counts_t)
AREG_exp_raw$AREG_RawCounts <- AREG_exp_raw$`total_counts_t$ENSG00000109321.11`
reorder_cells <- c("MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                   "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                   "T47D_rep1", "T47D_rep2", "T47D_rep3",
                   "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                   "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                   "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
                   "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                   "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                   "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
AREG_exp_raw <- AREG_exp_raw %>%
  dplyr :: slice(match(reorder_cells, cell_line))
AREG_exp_raw$cell_line <- factor(AREG_exp_raw$cell_line, levels = AREG_exp_raw$cell_line)
ggplot(AREG_exp_raw, aes(x = cell_line,
                     y = AREG_RawCounts)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_continuous(palette="magma")
  #scale_color_manual(values = c("red", "red", "red",
                               "blue", "blue", "blue",
                               "green", "green", "green",
                               "red", "red", "red",
                               "blue", "blue", "blue",
                               "green", "green", "green",
                               "red", "red", "red",
                               "blue", "blue", "blue",
                               "green", "green", "green"))

#IRS2
ETV5_exp <- as.data.frame(normalized_total_counts_t$ENSG00000185950.9)
ETV5_exp$cell_line <- rownames(normalized_total_counts_t)
ETV5_exp$ETV5_NormCounts <- ETV5_exp$`normalized_total_counts_t$ENSG00000185950.9`
reorder_cells <- c("MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                   "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                   "T47D_rep1", "T47D_rep2", "T47D_rep3",
                   "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                   "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                   "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
                   "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                   "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                   "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
ETV5_exp <- ETV5_exp %>%
  dplyr :: slice(match(reorder_cells, cell_line))
ETV5_exp$cell_line <- factor(ETV5_exp$cell_line, levels = ETV5_exp$cell_line)
ggplot(ETV5_exp, aes(x = cell_line,
                         y = ETV5_NormCounts)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))

#ETV5
ETV5_exp <- as.data.frame(normalized_total_counts_t$ENSG00000244405.8)
ETV5_exp$cell_line <- rownames(normalized_total_counts_t)
ETV5_exp$ETV5_NormCounts <- ETV5_exp$`normalized_total_counts_t$ENSG00000244405.8`
reorder_cells <- c("MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                   "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                   "T47D_rep1", "T47D_rep2", "T47D_rep3",
                   "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                   "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                   "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
                   "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                   "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                   "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
ETV5_exp <- ETV5_exp %>%
  dplyr :: slice(match(reorder_cells, cell_line))
ETV5_exp$cell_line <- factor(ETV5_exp$cell_line, levels = ETV5_exp$cell_line)
ggplot(ETV5_exp, aes(x = cell_line,
                     y = ETV5_NormCounts)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))

#THRB
THRB_exp <- as.data.frame(normalized_total_counts_t$ENSG00000151090.19)
THRB_exp$cell_line <- rownames(normalized_total_counts_t)
THRB_exp$THRB_NormCounts <- THRB_exp$`normalized_total_counts_t$ENSG00000151090.19`
reorder_cells <- c("MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                   "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                   "T47D_rep1", "T47D_rep2", "T47D_rep3",
                   "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                   "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                   "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
                   "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                   "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                   "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
THRB_exp <- THRB_exp %>%
  dplyr :: slice(match(reorder_cells, cell_line))
THRB_exp$cell_line <- factor(THRB_exp$cell_line, levels = THRB_exp$cell_line)
ggplot(THRB_exp, aes(x = cell_line,
                     y = THRB_NormCounts)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))

#PRLR
PRLR_exp <- as.data.frame(normalized_total_counts_t$ENSG00000113494.17)
PRLR_exp$cell_line <- rownames(normalized_total_counts_t)
PRLR_exp$PRLR_NormCounts <- PRLR_exp$`normalized_total_counts_t$ENSG00000113494.17`
reorder_cells <- c("MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                   "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                   "T47D_rep1", "T47D_rep2", "T47D_rep3",
                   "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                   "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                   "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
                   "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                   "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                   "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
PRLR_exp <- PRLR_exp %>%
  dplyr :: slice(match(reorder_cells, cell_line))
PRLR_exp$cell_line <- factor(PRLR_exp$cell_line, levels = PRLR_exp$cell_line)
ggplot(PRLR_exp, aes(x = cell_line,
                     y = PRLR_NormCounts)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90)) #+

#PIK3R1
PIK3R1_exp <- as.data.frame(normalized_total_counts_t$ENSG00000145675.15)
PIK3R1_exp$cell_line <- rownames(normalized_total_counts_t)
PIK3R1_exp$PIK3R1_NormCounts <- PIK3R1_exp$`normalized_total_counts_t$ENSG00000145675.15`
reorder_cells <- c("MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                   "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                   "T47D_rep1", "T47D_rep2", "T47D_rep3",
                   "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                   "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                   "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
                   "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                   "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                   "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
PIK3R1_exp <- PIK3R1_exp %>%
  dplyr :: slice(match(reorder_cells, cell_line))
PIK3R1_exp$cell_line <- factor(PIK3R1_exp$cell_line, levels = PIK3R1_exp$cell_line)
ggplot(PIK3R1_exp, aes(x = cell_line,
                     y = PIK3R1_NormCounts)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))

#PIK3R1
PIK3R1_exp <- as.data.frame(normalized_total_counts_t$ENSG00000145675.15)
PIK3R1_exp$cell_line <- rownames(normalized_total_counts_t)
PIK3R1_exp$PIK3R1_NormCounts <- PIK3R1_exp$`normalized_total_counts_t$ENSG00000145675.15`
reorder_cells <- c("MCF7_corrected_rep1", "MCF7_corrected_rep2", "MCF7_corrected_rep3",
                   "MCF7_parental_rep1", "MCF7_parental_rep2", "MCF7_parental_rep3",
                   "T47D_rep1", "T47D_rep2", "T47D_rep3",
                   "HTert_WT_rep1", "HTert_WT_rep2", "HTert_WT_rep3",
                   "HTert_ex9_rep1", "HTert_ex9_rep2", "HTert_ex9_rep3",
                   "HTert_ex20_rep1", "HTert_ex20_rep2", "HTert_ex20_rep3",
                   "MCF10A_parental_rep1", "MCF10A_parental_rep2", "MCF10A_parental_rep3",
                   "MCF10A_E545K_rep1", "MCF10A_E545K_rep2", "MCF10A_E545K_rep3",
                   "MCF10A_H1047R_rep1", "MCF10A_H1047R_rep2", "MCF10A_H1047R_rep3")
PIK3R1_exp <- PIK3R1_exp %>%
  dplyr :: slice(match(reorder_cells, cell_line))
PIK3R1_exp$cell_line <- factor(PIK3R1_exp$cell_line, levels = PIK3R1_exp$cell_line)
ggplot(PIK3R1_exp, aes(x = cell_line,
                       y = PIK3R1_NormCounts)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))





###Comparisons for doing shared pathway analyses
total_dds_cellSpec <- DESeqDataSetFromMatrix(countData= total_counts,
                                       colData= total_coldata_v2,
                                       design= ~ total_treatment_v2) 
total_dds_cellSpec <- DESeq(total_dds_cellSpec)
resultsNames(total_dds_cellSpec)
HTert_ex9vex20_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "ex9_HTert", "ex20_HTert"))
cancer_ex9vex20_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "ex9_MCF7", "ex20_T47D"))
MCF10A_ex9vex20_res <- results(total_dds_cellSpec, contrast = c("total_treatment_v2", "ex9_MCF10A", "ex20_MCF10A"))


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
dotplot(HTert_ex9vex20_KEGG,
        showCategory = 20,
        title = "HTert L2fc(+/-1.5) KEGG")
dotplot(cancer_ex9vex20_KEGG,
        showCategory = 20,
        title = "cancer L2fc(+/-1.5) KEGG")
dotplot(MCF10A_ex9vex20_KEGG,
        showCategory = 20,
        title = "MCF10A L2fc(+/-1.5) KEGG")

HTert_ex9vex20_KEGG_df <- as.data.frame(HTert_ex9vex20_KEGG)
cancer_ex9vex20_KEGG_df <- as.data.frame(cancer_ex9vex20_KEGG)
MCF10A_ex9vex20_KEGG_df <- as.data.frame(MCF10A_ex9vex20_KEGG)

sharedPaths <- intersect(intersect(HTert_ex9vex20_KEGG_df$ID,
                                   cancer_ex9vex20_KEGG_df$ID),
                         MCF10A_ex9vex20_KEGG_df$ID)
sharedPaths
sharedPaths_df <- HTert_ex9vex20_KEGG_df %>%
  subset(ID %in% sharedPaths) %>%
  select(ID, Description)
sharedPaths_df