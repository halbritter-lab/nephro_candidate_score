# Title: Positive kidney genes

# load libraries
library(tidyverse)
library(utils)

# download and unzip genes associated with kidney disease from github repository kidney-genetics
kg_url <- paste0("https://github.com/halbritter-lab/kidney-genetics/raw/main/analyses/A_MergeAnalysesSources/results/A_MergeAnalysesSources.", kidney_genetics_version, ".csv.gz")

download.file(url = kg_url,
              destfile = paste0("gene_score/labels/raw/A_MergeAnalysesSources.", kidney_genetics_version, ".csv.gz"))

gunzip(filename = paste0("gene_score/labels/raw/A_MergeAnalysesSources.", kidney_genetics_version, ".csv.gz"), 
       destname = paste0("gene_score/labels/raw/A_MergeAnalysesSources.", kidney_genetics_version, ".csv"))

# load data
kid_gen <- read.csv(paste0("gene_score/labels/raw/A_MergeAnalysesSources.", kidney_genetics_version, ".csv"))

# filter genes above a specific evidence count
evid_thresh <- 0 # TODO change
pos_genes <- kid_gen %>% 
  filter(evidence_count >= evid_thresh)


# TODO: selection of evidence group etc...
