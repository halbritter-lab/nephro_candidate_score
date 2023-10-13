# Title: Positive kidney genes

# load libraries
library(tidyverse)
library(R.utils)
library(config)

# read configs
config_vars <- config::get(file = "config.yml")
script_path <- "gene_score/labels"

# save current working directory
wd_bef_script_exe <- getwd()

# set new working directory
setwd(file.path(config_vars$PROJECT_DIR, script_path))

# download and unzip genes associated with kidney disease from github repository kidney-genetics
kg_url <- paste0("https://github.com/halbritter-lab/kidney-genetics/blob/main/analyses/A_MergeAnalysesSources/results/A_MergeAnalysesSources.", config_vars$kidney_genetics_version, ".csv.gz?raw=TRUE")

download.file(url = kg_url,
              destfile = paste0("raw/A_MergeAnalysesSources.", config_vars$kidney_genetics_version, ".csv.gz"))

gunzip(filename = paste0("raw/A_MergeAnalysesSources.", config_vars$kidney_genetics_version, ".csv.gz"), 
       destname = paste0("raw/A_MergeAnalysesSources.", config_vars$kidney_genetics_version, ".csv"))

# load data
kid_gen <- read.csv(paste0("raw/A_MergeAnalysesSources.", config_vars$kidney_genetics_version, ".csv")) %>% 
  rename(hgnc_id_int = "hgnc_id")


# write results
write.csv(kid_gen, paste0("results/positive_genes_", config_vars$creation_date, ".csv"), row.names = FALSE)

# set back former working directory
setwd(wd_bef_script_exe)
