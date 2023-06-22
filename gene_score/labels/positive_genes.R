# Title: Positive kidney genes

# load libraries
library(tidyverse)
library(utils)

# load genes associated with kidney disease from github repository kidney-genetics
kg_version <- "2023-05-18"  # TODO: change 
kg_url <- paste0("https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/merged/KidneyGenetics_MergeAnalysesSources.", kg_version, ".csv")
kid_gen <- read.csv(kg_url)

# filter genes above a specific evidence count
evid_thresh <- 0 # TODO change
pos_genes <- kid_gen %>% 
  filter(evidence_count >= evid_thresh)


# TODO: selection of evidence group etc...
