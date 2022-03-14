#DRUMLR
# devtools::install_github(repo="CutillasLab/DRUMLR", subdir="DRUMLR_code")
library(DRUMLR)
library(tidyverse)
if (!requireNamespace("BiocManager", quietly=TRUE))
# BiocManager::install('PharmacoGx')
library(PharmacoGx)

#check available datasets in pharmacoGx
gxsets <- availablePSets()
CCLE_Pset <- downloadPSet("CCLE_2015")
GDSC_Pset <- downloadPSet("GDSC_2019(v2-8.0)")

#build model using pharmacoGx pset
Model_info <- DRUMLR::BuildDRUMLs(CCLE_Pset, df_response = df_response)


DRUMLR::