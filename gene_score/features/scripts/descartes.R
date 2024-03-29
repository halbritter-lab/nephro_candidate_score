# Title: Descartes - ssRNA fetal kidney

# load libraries
library(tidyverse)	
library(jsonlite) 
library(config)
library(R.utils)

# define relative script path
project_topic <- "nephrology"
project_name <- "nephro_candidate_score"
script_path <- "/gene_score/features/"

# read configs
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"),
                           config = project_topic)

# set working directory
setwd(paste0(config_vars$projectsdir, project_name, script_path))

# download and read in main clusters and gene IDs
main_clusters_url <- "https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/labels/kidney_main_cluster_name_lbls.json"

download.file(main_clusters_url, 
              destfile = paste0("raw/descartes_fetal_kidney_main_clusters_", config_vars$creation_date, ".json"))

main_clusters <- fromJSON( paste0("raw/descartes_fetal_kidney_main_clusters_", config_vars$creation_date, ".json")) 

gene_ids_url <- "https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/counts/gene2ens.csv"
download.file(gene_ids_url, 
              destfile = paste0("raw/descartes_gene_ids_url_", config_vars$creation_date, ".csv"))

gene_ids <- read.csv(paste0("raw/descartes_gene_ids_url_", config_vars$creation_date, ".csv"), header=FALSE) %>% 
  dplyr::rename(symbol = V1, ensembl_gene_id = V2)


# function for downloading TPM values and cell percentage for each kidney cell type
get_cell_tpm <- function(cell_type){
  download_url <- paste0("https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/tables/cell_tpm/kidney/", cell_type, ".csv")
  download.file(download_url, 
                destfile = paste0("raw/descartes_", cell_type, "_tpm_", config_vars$creation_date, ".csv"))
}

get_cell_percentage <- function(cell_type){
  download_url <- paste0("https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/tables/cell_percentage/kidney/", cell_type, ".csv")
  download.file(download_url, 
                destfile = paste0("raw/descartes_", cell_type, "_percentage_", config_vars$creation_date, ".csv"))
}

# kidney cell types
cell_types <- c("mesangial", "metanephric", "ureteric_bud", "stromal", "vascular_endothelial")

# create dataframes for TPM values and cell expression data
perc_expr_df <- gene_ids %>% 
  filter(ensembl_gene_id %in% HGNC_table$ensembl_gene_id)

tpm_df <- gene_ids %>% 
  filter(ensembl_gene_id %in% HGNC_table$ensembl_gene_id)

# fill dataframes
for (cell_type in cell_types){
  get_cell_percentage(cell_type)
  get_cell_tpm(cell_type)
  
  ## percent expression
  ct_perc <- read.csv(paste0("raw/descartes_", cell_type, "_percentage_", config_vars$creation_date, ".csv"), header=FALSE)

  # in case multiple values are available for one gene, take the highest number
  ct_perc <- ct_perc %>% 
    group_by(V1) %>% 
    summarise(V2 = max(V2, na.rm = TRUE))
  colnames(ct_perc) <- c("symbol", paste0(cell_type, "_perc_expr"))
  
  perc_expr_df <- left_join(perc_expr_df, ct_perc, by = "symbol")
  
  ## tpm values
  ct_tpm <- read.csv(paste0("raw/descartes_", cell_type, "_tpm_", config_vars$creation_date, ".csv"), header=FALSE)
  
  # in case multiple values are available for one gene, take the highest number
  ct_tpm <- ct_tpm %>% 
    group_by(V1) %>% 
    summarise(V2 = max(V2, na.rm = TRUE))
  colnames(ct_tpm) <- c("symbol", paste0(cell_type, "_tpm"))
  
  tpm_df <- left_join(tpm_df, ct_tpm, by = "symbol")
}

# perc_expr_df <- perc_expr_df %>% dplyr::rename(symbol = symbol)
# tpm_df <- tpm_df %>% dplyr::rename(symbol = symbol)

# write results
write.csv(perc_expr_df, 
          paste0("results/descartes_fetal_kidney_percent_expression_", config_vars$creation_date, ".csv"), 
          row.names = FALSE)

gzip(paste0("results/descartes_fetal_kidney_percent_expression_", config_vars$creation_date, ".csv"),
     overwrite = TRUE)

write.csv(tpm_df, 
          paste0("results/descartes_fetal_kidney_tpm_", config_vars$creation_date, ".csv"), 
          row.names = FALSE)

gzip(paste0("results/descartes_fetal_kidney_tpm_", config_vars$creation_date, ".csv"),
     overwrite = TRUE)

## Calculate nTPM values
# Calculate the effective library sizes on complete rows
tpm_df_cc <- tpm_df %>% column_to_rownames(var = "ensembl_gene_id") %>% dplyr::select(-symbol)
tpm_df_cc <- tpm_df_cc[complete.cases(tpm_df_cc),] 

# function for quantile normalization (https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/)
quantile_normalisation <- function(df){
  df_rank <- apply(df, 2, rank, ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

# perform quantile normalization on TPM complete cases dataframe
ntpm_df_cc <- quantile_normalisation(tpm_df_cc)
colnames(ntpm_df_cc) <- gsub("_tpm", "_nTPM", colnames(ntpm_df_cc))
ntpm_df_cc_with_id <- ntpm_df_cc %>% data.frame() %>% rownames_to_column(var = "ensembl_gene_id")

# calculate the normalized tissue specifitiy index tau according to Yanai et al. 
N <- ncol(ntpm_df_cc)
max_row <- apply(ntpm_df_cc, 1, max)
norm_data <- ntpm_df_cc / max_row
tau_values <- sapply(1:nrow(norm_data), function(i) {
  xi <- norm_data[i, ]
  sum((1 - xi) / (N - 1))
})

fetal_kidney_tau_df <- data.frame(ensembl_gene_id = rownames(norm_data), fetal_kidney_tau = tau_values) 

# add symbol column
fetal_kidney_tau_df <- fetal_kidney_tau_df %>% left_join(distinct(tpm_df[, c("ensembl_gene_id", "symbol")]), by = "ensembl_gene_id")
ntpm_df_cc_with_id <- ntpm_df_cc_with_id %>% left_join(distinct(tpm_df[, c("ensembl_gene_id", "symbol")]), by = "ensembl_gene_id")

# write results
write.csv(fetal_kidney_tau_df, 
          paste0("results/descartes_fetal_kidney_tau_", config_vars$creation_date, ".csv"), 
          row.names = FALSE)

gzip(paste0("results/descartes_fetal_kidney_tau_", config_vars$creation_date, ".csv"),
     overwrite = TRUE)

write.csv(ntpm_df_cc_with_id, 
          paste0("results/descartes_fetal_nptm_", config_vars$creation_date, ".csv"), 
          row.names = FALSE)

gzip(paste0("results/descartes_fetal_nptm_", config_vars$creation_date, ".csv"), 
     overwrite = TRUE)

