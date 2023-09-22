# Title: Descartes - ssRNA fetal kidney

#https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/mesangial/in/kidney



# load libraries
library(tidyverse)	
library(jsonlite) 


# download and read in files
main_clusters_url <- "https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/labels/kidney_main_cluster_name_lbls.json"

download.file(main_clusters_url, 
              destfile = paste0("gene_score/features/raw/descartes_fetal_kidney_main_clusters_", creation_date, ".json"))

main_clusters <- fromJSON( paste0("gene_score/features/raw/descartes_fetal_kidney_main_clusters_", creation_date, ".json")) 

gene_ids_url <- "https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/counts/gene2ens.csv"
download.file(gene_ids_url, 
              destfile = paste0("gene_score/features/raw/descartes_gene_ids_url_", creation_date, ".csv"))

gene_ids<- read.csv(paste0("gene_score/features/raw/descartes_gene_ids_url_", creation_date, ".csv"), header=FALSE) %>% 
  rename(gene_symbol = V1, ensembl_gene_id = V2)


# download tpm and cell percentage for each kidney cell type
get_cell_tpm <- function(cell_type){
  download_url <- paste0("https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/tables/cell_tpm/kidney/", cell_type, ".csv")
  download.file(download_url, 
                destfile = paste0("gene_score/features/raw/descartes_", cell_type, "_tpm_", creation_date, ".csv"))
}

get_cell_percentage <- function(cell_type){
  download_url <- paste0("https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/tables/cell_percentage/kidney/", cell_type, ".csv")
  download.file(download_url, 
                destfile = paste0("gene_score/features/raw/descartes_", cell_type, "_percentage_", creation_date, ".csv"))
}

# get_cell_tpm("mesangial")
# get_cell_percentage("mesangial")
# get_cell_tpm("metanephric")
# get_cell_percentage("metanephric")
# get_cell_tpm("uretic_bud")
# get_cell_percentage("uretic_bud")
# get_cell_tpm("stromal")
# get_cell_percentage("stromal")
# get_cell_tpm("vascular_endothelial")
# get_cell_percentage("vascular_endothelial")


cell_types <- c("mesangial", "metanephric", "ureteric_bud", "stromal", "vascular_endothelial")

perc_expr_df <- gene_ids %>% 
  filter(ensembl_gene_id %in% HGNC_table$ensembl_gene_id)

tpm_df <- gene_ids %>% 
  filter(ensembl_gene_id %in% HGNC_table$ensembl_gene_id)

for (cell_type in cell_types){
  get_cell_percentage(cell_type)
  get_cell_tpm(cell_type)
  
  ## percent expression
  ct_perc <- read.csv(paste0("gene_score/features/raw/descartes_", cell_type, "_percentage_", creation_date, ".csv"), header=FALSE)

  # in case multiple values are available for one gene, take the highest number
  ct_perc <- ct_perc %>% 
    group_by(V1) %>% 
    summarise(V2 = max(V2, na.rm = TRUE))
  colnames(ct_perc) <- c("gene_symbol", paste0(cell_type, "_perc_expr"))
  
  perc_expr_df <- left_join(perc_expr_df, ct_perc, by = "gene_symbol")
  
  ## tpm values
  ct_tpm <- read.csv(paste0("gene_score/features/raw/descartes_", cell_type, "_tpm_", creation_date, ".csv"), header=FALSE)
  
  # in case multiple values are available for one gene, take the highest number
  ct_tpm <- ct_tpm %>% 
    group_by(V1) %>% 
    summarise(V2 = max(V2, na.rm = TRUE))
  colnames(ct_tpm) <- c("gene_symbol", paste0(cell_type, "_tpm"))
  
  tpm_df <- left_join(tpm_df, ct_tpm, by = "gene_symbol")
}

perc_expr_df <- perc_expr_df %>% select(-gene_symbol)
tpm_df <- tpm_df %>% select(-gene_symbol)

# write results
write.csv(perc_expr_df, paste0("gene_score/features/results/descartes_fetal_kidney_percent_expression_", creation_date, ".csv"), row.names = FALSE)
write.csv(tpm_df, paste0("gene_score/features/results/descartes_fetal_kidney_tpm_", creation_date, ".csv"), row.names = FALSE)





