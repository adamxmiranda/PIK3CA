library(tidyverse)
library(readr)
library(pheatmap)
library(gplots)
library(ggplot2)
library(factoextra)
library(cluster)

setwd("~/Desktop/PhD_Research/WGBS_MethPipe_Islet_Cells_Copy/WGBS_R/K-means Heatmap")

##use paste command at the command line to merge meth columns in parallel from each individual into one master file## 
##this matrix should contain the following columns: chr, start, end, 1A, 2A, 3B, 4B## 

# load in matrix listing avg methyl scores of HMRs for individuals 
# methyl scores were calculated using bedtools map OR roimethstat function (MethPipe) using the methcounts file that only lists methylation levels for symmetric and mutated CpG sites
# avg methyl scores were calculated with bedtools map
avg_hmr_matrix_CpG <- read_tsv("matrix_avg_meth_hmr_CpG.bed", col_names = c("Chr","Start","End","1A","2A","3B","4B"))
avg_hmr_matrix_CpG_cleaned <- avg_hmr_matrix_CpG %>%
  dplyr::select("1A","2A","3B","4B") 

avg_hmr_matrix_CpG_roimethstat <- read_tsv("matrix_avg_meth_hmr_CpG.bed", col_names = c("Chr","Start","End","1A","2A","3B","4B")) 
avg_hmr_matrix_CpG_roimethstat_cleaned <- avg_hmr_matrix_CpG_roimethstat %>%
  dplyr::select("1A","2A","3B","4B") 

# visualize how hmr methyl levels are distributed in the dataset usig geom_histogram() and geom_density()
# methyl level matrix created from roimethstat and bedtools map have the same output
ggplot(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG_roimethstat$`1A`)) +
  geom_histogram(col=I("red"), alpha=.8) + 
  geom_histogram(data = avg_hmr_matrix_CpG, aes(avg_hmr_matrix_CpG$`1A`), col=I("blue"), alpha=.5)

ggplot(data = avg_hmr_matrix_CpG, aes(avg_hmr_matrix_CpG$`1A`)) +
  geom_histogram(col=I("blue"), alpha=.8)

# density curve
# visualizes the underlying probability distribution of the data by drawing an appropriate continuous curve 
ggplot(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG_roimethstat$`1A`)) +
  geom_density()

plot <- ggplot(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG_roimethstat$`1A`)) +
  geom_density(color="darkblue") + 
  geom_density(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG$`2A`), color="darkgreen") +
  geom_density(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG$`3B`), color="darkred") +
  geom_density(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG$`4B`))
plot + labs(title = "Average Methylation Level Density Curve for HMRs ")


d1 <- ggplot(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG_roimethstat$`1A`)) +
  geom_density(color="darkblue")
d1

d2 <- ggplot(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG_roimethstat$`2A`)) +
  geom_density(color="darkgreen")
d2

d3 <- ggplot(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG_roimethstat$`3B`)) +
  geom_density(color="darkred")
d3

d4 <- ggplot(data=avg_hmr_matrix_CpG_roimethstat, aes(avg_hmr_matrix_CpG_roimethstat$`4B`)) +
  geom_density(color="purple")
d4


# Does any member of the variable x have the value NA?”
# 1A, 2A, and 4B have have NA valeus 
apply(avg_hmr_matrix_CpG_cleaned, 2, function(x) any(is.na(x)))
# Does any member of the variable x have the value NaN (Not a Number)?”
# all false 
apply(avg_hmr_matrix_CpG_cleaned, 2, function(x) any(is.nan(x)))
# Does any member of the variable x have the value -Inf or Inf?”
# all false 
apply(avg_hmr_matrix_CpG_cleaned, 2, function(x) any(is.infinite(x)))
# Does any member of the variable x have the value NA or -Inf or Inf?”
# applies both functions 
apply(avg_hmr_matrix_CpG_cleaned, 2, function(x) any(is.na(x) | is.infinite(x)))

# count missing values
# number of missing values is 3 
missing_values_avg_hmr_matrix_CpG_cleaned <- avg_hmr_matrix_CpG_cleaned %>%
  dplyr::summarise(count = sum(is.na(avg_hmr_matrix_CpG_cleaned)))

# count missing values and distinct values (which also counts missing values)
df <- avg_hmr_matrix_CpG_cleaned %>%
  dplyr::summarise(n = n_distinct(avg_hmr_matrix_CpG_cleaned), 
                   na = sum(is.na(avg_hmr_matrix_CpG_cleaned))) 

# converts data into a tibble (i.e., tbl_df)
# my_data <- as_tibble(avg_hmr_matrix_CpG_cleaned)
# drop_na(my_data)

# drop rows containg missing values
# why doesn't the function below work? 
avg_hmr_matrix_CpG_cleaned_v2 <- avg_hmr_matrix_CpG_cleaned %>%
  tidyr::drop_na(avg_hmr_matrix_CpG_cleaned)

# drop rows containg missing values
# only 1 row is dropped
drop_na(avg_hmr_matrix_CpG_cleaned)
avg_hmr_matrix_CpG_cleaned_v2 <- avg_hmr_matrix_CpG_cleaned %>% drop_na()

# load in listing avg methyl scores for individuals using the methcounts file that lists methylation for all C's, including CpG C's and others 
# avg methyl scores were calculated with bedtools map 
avg_hmr_matrix_all <- read_tsv("matrix_avg_meth_hmr_all.bed", col_names = c("Chr","Start","End","1A","2A","3B","4B"))
avg_hmr_matrix_all_cleaned <- avg_hmr_matrix_all %>%
  dplyr::select("1A","2A","3B","4B")

# Does any member of the variable x have the value NA (Not Available)?”
# is.na detects both values: Na and NaN 
# 1A and 2A have NA values 
apply(avg_hmr_matrix_all_cleaned, 2, function(x) any(is.na(x)))
# Does any member of the variable x have the value NaN (Not a Number)?”
# all false 
apply(avg_hmr_matrix_all_cleaned, 2, function(x) any(is.nan(x)))
# Does any member of the variable x have the value -Inf or Inf?”
# all false 
apply(avg_hmr_matrix_all_cleaned, 2, function(x) any(is.infinite(x)))
# Does any member of the variable x have finite values (i.e., not infinite and not missing)?”
# all true 
apply(avg_hmr_matrix_all_cleaned, 2, function(x) any(is.finite(x)))
# Does any member of the variable x have the value NA or -Inf or Inf?”
# applies both functions 
apply(avg_hmr_matrix_all_cleaned, 2, function(x) any(is.na(x) | is.infinite(x)))

# counting missing values
# number of missing values is 2
missing_values_avg_hmr_matrix_all_cleaned <- avg_hmr_matrix_all_cleaned %>%
  dplyr::summarise(count = sum(is.na(avg_hmr_matrix_all_cleaned)))

# count missing values and distinct values (which also counts missing values)
df2 <- avg_hmr_matrix_all_cleaned%>%
  dplyr::summarise(n = n_distinct(avg_hmr_matrix_all_cleaned), 
                   na = sum(is.na(avg_hmr_matrix_all_cleaned))) 

# drop rows containg missing values
drop_na(avg_hmr_matrix_all_cleaned)
avg_hmr_matrix_all_cleaned_v2 <- avg_hmr_matrix_all_cleaned %>% drop_na()

#indx <- apply(avg_hmr_matrix_CpG_cleaned, 2, function(x) any(is.na(x) | is.infinite(x)))
#colnames[indx]

# pheatmap function 
# default behavior of the function includes the hierarchical clustering of both rows and columns
# objects within clusters are as similar as possible (i.e., high intra-class similarity)
# in k-means clustering, each cluster is represented by its center (i.e., centorid) which corresponds to the mean of points assigned to the cluster
# scale is Eucledian distance 
# if you want to turn off the clustering, you can set either cluster_cols or cluster_rows to False
pheatmap(avg_hmr_matrix_all_cleaned_v2, kmeans_k = 4)
p <- pheatmap(normalized_total_counts, kmeans_k = 9)
clusters <- as.data.frame(p$kmeans$cluster)
# distance measure used in clustering rows
# possible values are "correlation" for Pearson correlation and all the distances supported by dist, such as "euclidean", etc.
# if the value is none of the above it is assumed that a distance matrix is provided
pheatmap(avg_hmr_matrix_all_cleaned_v2, clustering_distance_rows= "correlation", clustering_distance_cols = "correlation")
# by setting kmeans_k to NA, rows are not aggregated 
pheatmap(avg_hmr_matrix_all_cleaned_v2, kmeans_k = NA, clustering_method = "ward.D2")

pheatmap(avg_hmr_matrix_CpG_cleaned_v2, kmeans_k = 4) 
pheatmap(avg_hmr_matrix_Cpg_cleaned_v2, kmeans_k = NA, clustering_method = "ward.D2")
pheatmap(normalized_total_counts, clustering_distance_rows= "correlation", clustering_distance_cols = "correlation")

# compute gap statistic 
gap_stat <- clusGap(avg_hmr_matrix_CpG_cleaned_v2, FUN = kmeans, nstart = 30, K.max = 10, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")

#elbow method
kmax <- 15
wss <-sapply(1:kmax,
             function(k){kmeans(normalized_total_counts, k,
                                nstart = 50, iter.max= 15) $tot.withinss})
wss
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# identify optimal k-means value
# how to choose the right number of expected clusters (k)?
# the idea is to compute k-means clustering using different values of clusters k
kmeans(avg_hmr_matrix_all_cleaned_v2, 4, iter.max = 10, nstart = 1)

# average silhouette for kmeans 
# determines and visualizes the optimal number of clusters using silhouette 
fviz_nbclust(avg_hmr_matrix_all_cleaned_v2, kmeans, method = "silhouette") 

# if you can't run locally due to not enough RAM 
mclapply(avg_hmr_matrix_all_cleaned_v2, fviz_nbclust) 

# gap statistic
library(cluster)
# compute gap statistic for k-means
# recommended valie is ~500
gap_stat <- clusGap(iris.scaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 10)
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)

# heatmap.2
mat <- matrix(avg_hmr_matrix_all_cleaned_v2)
heatmap.2(mat)

