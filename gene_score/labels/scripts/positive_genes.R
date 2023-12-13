# Title: Positive kidney genes

# load libraries
library(tidyverse)
library(R.utils)
library(config)

# define relative script path
project_topic <- "nephrology"
project_name <- "nephro_candidate_score"
script_path <- "/gene_score/labels/"

# read configs
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"),
                           config = project_topic)

# set working directory
setwd(paste0(config_vars$projectsdir, project_name, script_path))

# download genes associated with kidney disease from github repository "kidney-genetics"
destfile <-  paste0("raw/A_MergeAnalysesSources.", config_vars$kidney_genetics_version, ".csv.gz")

if (!file.exists(destfile)) {
  # if the file doesn't exist, download it
  kg_url <- paste0("https://github.com/halbritter-lab/kidney-genetics/blob/main/analyses/A_MergeAnalysesSources/results/A_MergeAnalysesSources.", config_vars$kidney_genetics_version, ".csv.gz?raw=TRUE")
  download.file(url = kg_url,
                destfile = destfile)
  } else {
    cat("'", strsplit(destfile, "/")[[1]][-1], "' already exists. No need to download.", sep ="")
}

# load data
kid_gen <- read_csv(destfile, show_col_types = FALSE) %>% 
  rename(hgnc_id_int = hgnc_id)

# write results
write.csv(kid_gen, paste0("results/positive_genes_", config_vars$creation_date, ".csv"), 
          row.names = FALSE)

gzip(paste0("results/positive_genes_", config_vars$creation_date, ".csv"),
     overwrite = TRUE)

