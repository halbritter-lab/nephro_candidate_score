# Title: GTEx expression values of different tissues

# load libraries
library(tidyverse)
library(R.utils)
library(config)

# read configs
config_vars <- config::get(file = "config.yml")
script_path <- "gene_score/features"

# save current working directory
wd_bef_script_exe <- getwd()

# set working directory
setwd(file.path(config_vars$PROJECT_DIR, script_path))

# download and unzip GTEx RNA expression data of different tissues
gtex_download_url <- "https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip"
download.file(gtex_download_url,
              destfile = paste0("raw/rna_tissue_gtex_", config_vars$creation_date, ".tsv.zip"))

unzip(zipfile = paste0("raw/rna_tissue_gtex_", config_vars$creation_date, ".tsv.zip"),
      exdir = "raw/")  

unzip(zipfile = paste0("raw/rna_tissue_gtex_", creation_date, ".tsv.zip"),
      exdir = "raw/")  

file.rename("raw/rna_tissue_gtex.tsv", paste0("raw/rna_tissue_gtex_", config_vars$creation_date, ".tsv"))

# load data
rna_tissue_gtex_nTPM <- read.delim(paste0("raw/rna_tissue_gtex_", config_vars$creation_date, ".tsv")) %>% 
  dplyr::select(ensembl_gene_id = Gene, Tissue, nTPM) %>% 
  spread(key = Tissue, value = nTPM) 

symbol_df <- read.delim(paste0("raw/rna_tissue_gtex_", config_vars$creation_date, ".tsv")) %>% 
  dplyr::select(ensembl_gene_id = Gene, symbol = Gene.name) %>% 
  distinct()

# dataset contains expression values for all tissues for 19764 genes - except for retina: here for 20090 genes
# => remove genes, that only have expression values for retina, but no other tissues
rna_tissue_gtex_nTPM <- rna_tissue_gtex_nTPM[complete.cases(rna_tissue_gtex_nTPM), ]

# aggregate brain tissues
brain_regions <- c("amygdala", "caudate", "cerebellum", "cerebral cortex", "hippocampus", "hypothalamus", "nucleus accumbens", "putamen", "substantia nigra")
brain_nTPM <- rna_tissue_gtex_nTPM %>%
  dplyr::select(ensembl_gene_id, all_of(brain_regions))
brain_nTPM_med <- data.frame(ensembl_gene_id = brain_nTPM$ensembl_gene_id, brain_median = apply(brain_nTPM[, -1], 1, median))

# join median brain nTPM values with other tissue df
rna_tissue_gtex_nTPM_agg <- rna_tissue_gtex_nTPM %>% 
  dplyr::select(-all_of(brain_regions)) %>% 
  left_join(brain_nTPM_med, by = "ensembl_gene_id") %>%
  column_to_rownames(var = "ensembl_gene_id")

# calculate the normalized tissue specifitiy index tau according to Yanai et al. for all genes
N <- ncol(rna_tissue_gtex_nTPM_agg)
max_row <- apply(rna_tissue_gtex_nTPM_agg, 1, max)
norm_data <- rna_tissue_gtex_nTPM_agg / max_row
tau_values <- sapply(1:nrow(norm_data), function(i) {
  xi <- norm_data[i, ]
  sum((1 - xi) / (N - 1))
})

tau_df <- data.frame(ensembl_gene_id = rownames(norm_data), gtex_tau = tau_values) %>% 
  left_join(symbol_df, by = "ensembl_gene_id")

# add prefix and get rownames back to column
names(rna_tissue_gtex_nTPM_agg) <- paste0("gtex_", names(rna_tissue_gtex_nTPM_agg))

rna_tissue_gtex_nTPM_agg <- rna_tissue_gtex_nTPM_agg %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  rename_all(~ str_replace_all(.x, " ", "_")) %>% 
  left_join(symbol_df, by = "ensembl_gene_id")

# write results - nTPM values
write.csv(rna_tissue_gtex_nTPM_agg, paste0("results/rna_tissue_gtex_nTPM_agg_" , config_vars$creation_date, ".csv"), row.names = FALSE)

# write results - tau values
write.csv(tau_df, paste0("results/rna_tissues_gtex_nTPM_agg_tau_val_" , config_vars$creation_date, ".csv"), row.names = FALSE)

# set back former working directory
setwd(wd_bef_script_exe)

