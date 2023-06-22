# Title: Positive kidney genes

# load libraries
library(tidyverse)
library(utils)

# load genes associated with kidney disease from github repository kidney-genetics
kg_url <– paste0("https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/merged/KidneyGenetics_MergeAnalysesSources.", kidney-genetics_version, ".csv")
kid_gen <– read.csv(kg_url)
