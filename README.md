# Functional Genomic Analysis of *PIK3CA*  Multi-lineage Isogenic Breast Cell Model

This repo contains my code for functional genomic analysis of my multi-lineage *PIK3CA* isogenic cell line model.

## For analyses specifically included in the publication: ID# *pending* please see the folder titled Publication



## The rest of this repo contains code from the project at large. The preprocessing of fastq files is used in the analyses used in the publication. Some data/cell lines are not included in the final analyses for the paper.


## Cell Line Model

For this inquiry, we make use of a multi-lineage isogenic cell line model. In this model, we have cell lines of three different breast cell lineages representing each of the three different PIK3CA genotypes of interest(WT, E545K, H1047R). A breakdown of the lineages and genotypes used in our study can be found in the table below:


|Genotype &rarr; <br /> Lineage &darr;  | Wild Type | E545K  | H1047R |
|:----|  :----:  |:----: | :----:  |
|MCF10A| **MCF10A** Parental | **MCF10A** AAV modified | **MCF10A** AAV modified|
|H-Tert IMEC| **H-Tert** Parental | **H-Tert** AAV modified | **H-Tert** AAV modified|
|Cancer Cells| **MCF7** AAV modified | **MCF7** Parental | **T47D**|

## ATAC-Seq

We performed ATAC-seq for two replicates of each of the cell lines in our model. These samples were prepared according to the Hodges Lab ATAC-seq protocol and were sequeced by VANTAGE at ~100 million reads per sample.

### Trimming
As is the case with other sequencing methods, adapter sequences need to be trimmed from output reads. Our ATAC protocol relies on a Nextera based adapter system. In our pipeline we utilize trim-galore, a wrapper for cutadapt and fastqc, which searches for known adapter sequences, including Nextera, and removes them from reads.

[Trim Galore](https://github.com/FelixKrueger/TrimGalore)

```bash
ATAC_Trim_loop_1.slrm
ATAC2_Trim_loop_1.slrm
```

### Mapping
singleton reads are removed and reads are sorted using the repair.sh function in the bbmap package.
[BBMap](https://sourceforge.net/projects/bbmap/)

Alignment of reads to the hg38 genome was performed using the Burrows-Wheeler Aligner.
[BWA](http://bio-bwa.sourceforge.net/)

Reads are then filtered (map quality of 40, mitochondrial DNA, and Blacklist Regions)

```bash
ATAC_BWA_MapnClean_loop_2.slrm
ATAC2_BWA_MapnClean_loop_2.slrm
```

### Remove PCR Duplicates
PCR duplicates are removed from the bam files using PICARD
[Picard](https://broadinstitute.github.io/picard/)

```bash
ATAC_remove_dupes_3.slrm
ATAC2_remove_dupes_3.slrm
```
### Call Peaks

As is the case with other ATAC based methods, we identify accumulation of reads, and thus accessible regions, using peak calling methods. The method we use is Genrich, as recommended by Harvard FAS Informatics. Genrich incorporates all replicates initially into its peak calling algorithm. We prefer to use Genrich as it includes an ATAC read shift correction and its handling of biological replicates is more streamlined. However, different methods may be appropriate depending on your individual study.

[Genrich](https://github.com/jsh58/Genrich)

```bash
ATAC_Genrich_6_LG.slrm
ATAC2_Genrich_6_LG.slrm
```

### Further analyses
Further analyses can be found within my jupyter notebook
```bash
ATAC-seq.ipynb
```
## RNA-seq

We performed RNA-seq for three replicates of each of the cell lines in our model. RNA was prepared using the Qiagen RNeasy kit. Libraries were prepared and sequenced by VANTAGE using their ribo-depletion protocol at ~50 million reads per sample.

### Trimming and mapping

Our trimming step was performed as above using Trim Galore.

Mapping to the hg38 genome was performed using STAR Aligner
[STAR](https://github.com/alexdobin/STAR)

```bash
6142_RNA_seq_preprocess.slrm
```

##Filtering

Filtering of mapped reads was performed using SAMtools for a map quality of 30

##Counting

Reads were counted to genes using the featureCounts utility of the SubRead package. Gencode v32 gene coordinates were used.

[Subread](http://subread.sourceforge.net/)

```bash
6142_RNA_seq_featureCounts.slrm
```

###Further Analyses can be found in the following files

```bash
all_kmeans_heatmap.R
ClusterProfiler_MCFTrio.R
Total_DESeq.R
DRUML.R
GSVA_enrichment.R
LFC_Clustering.R
Lineage_Clustering.R
MCF7_HiC_overlap.R
SharedPaths.R
Total_clusterProfiler.R
Under_Over.R

###These are probably the most useful
Total_RNA_Pathway_Specific.ipynb
```
##Further Analyses
Analyses looking to combine the results of the RNA-seq and ATAC-seq data can be found in

```bash
DATA_Unite.ipynb
```

## Hi-C

Bead capture Hi-C was also performed on cells in our model. Many cell lines have yet to be prepared and this section is very much *under construction*.

These samples were prepared using a custom protocol and sequenced by VANTAGE at ~150 million reads per sample. Protocols.io link coming soon!
