# Title: Gene score data preprocessing master script

# load libraries
library(tidyverse)
library(utils)
library(R.utils)
library(config)

# read configs
config_vars <- config::get(file = "config.yml")
script_path <- "gene_score"

# set working directory
setwd(file.path(config_vars$PROJECT_DIR, script_path))

# source additional functions
source("hgnc_functions.R")
source("helper_functions.R")

# download HGNC gene table from github repository "kidney-genetics"
gene_table_url <- paste0("https://github.com/halbritter-lab/kidney-genetics/blob/main/analyses/B_AnnotationHGNC/results/non_alt_loci_set_coordinates.", config_vars$hgnc_gt_version, ".csv.gz?raw=TRUE")
download.file(url = gene_table_url, destfile = paste0("raw/HGNC_", config_vars$hgnc_gt_version, ".csv.gz"))
gunzip(filename = paste0("raw/HGNC_", config_vars$hgnc_gt_version, ".csv.gz"), 
       destname = paste0("raw/HGNC_", config_vars$hgnc_gt_version, ".csv"))

# load HGNC gene table and filter for protein-coding genes # TODO: change paths
HGNC_table <- read.csv(paste0("raw/HGNC_", config_vars$hgnc_gt_version, ".csv")) %>% 
  filter(locus_group == "protein-coding gene") %>% 
  dplyr::select(hgnc_id, entrez_id, ensembl_gene_id, symbol, alias_symbol, prev_symbol) %>% 
  rowwise() %>% 
  mutate(hgnc_id_int = as.integer(str_replace(hgnc_id, "HGNC:", "")))

# extract all protein coding gene symbols
all_prot_coding_gene_symbols <- unlist(strsplit(c(HGNC_table$symbol, HGNC_table$alias_symbol, HGNC_table$prev_symbol), "\\|")) %>% unique()


##### LABELS ##### 
# get positive genes
cat("get positive genes...")
source("labels/scripts/positive_genes.R")

# get dipensible genes #TODO: test with stable internet connection
cat("get dispensible genes...")
source("labels/scripts/dispensable_genes.R")


##### FEATURES ##### 
# get cellXgene features and join with HGNC table
cat("get cellXgene features...")
source("features/scripts/cellxgene.R")
cellXgene_ds1 <- read_csv(paste0("gene_score/features/results/cellxgene_expr_0b4a15a7-4e9e-4555-9733-2423e5c66469_", creation_date, ".csv"), show_col_types = FALSE)
# cellXgene_ds2 <- read_csv(paste0("gene_score/features/results/cellxgene_expr_d7dcfd8f-2ee7-4385-b9ac-e074c23ed190_", creation_date, ".csv"), show_col_types = FALSE)

HGNC_table <- HGNC_table %>%
  left_join(cellXgene_ds1, by = "ensembl_gene_id") #%>%
  # left_join(cellXgene_ds2, by = "ensembl_gene_id")


# get Descartes features and join with HGNC table
cat("get descartes fetal kidney features...")
source("gene_score/features/scripts/descartes.R")
descartes_fetal_kid_tau <- read_csv(paste0("gene_score/features/results/descartes_fetal_kidney_tau_", creation_date, ".csv"), show_col_types = FALSE)
descartes_fetal_kid_pe <- read_csv(paste0("gene_score/features/results/descartes_fetal_kidney_percent_expression_", creation_date, ".csv"), show_col_types = FALSE)
descartes_fetal_kid_nptm <- read_csv(paste0("gene_score/features/results/descartes_fetal_nptm_", creation_date, ".csv"), show_col_types = FALSE)

HGNC_table <- HGNC_table %>%
  left_join(descartes_fetal_kid_tau, by = "ensembl_gene_id") %>% 
  left_join(descartes_fetal_kid_pe, by = "ensembl_gene_id") %>% 
  left_join(descartes_fetal_kid_nptm, by = "ensembl_gene_id")
  

# get gnomAD features and join with HGNC table
cat("get gnomad features...")
source("gene_score/features/scripts/gnomad.R")
gnomad_constraints <- read_csv(paste0("gene_score/features/results/gnomad_constraints_", creation_date, ".csv"), show_col_types = FALSE)
# HGNC_table <- HGNC_table %>% left_join(gnomad_constraints, by = "ensembl_gene_id")
# HGNC_table <- left_join_rescue(HGNC_table, gnomad_constraints, by1 = "ensembl_gene_id", by2 = "symbol")
HGNC_table <- left_join_rescue_symbol(HGNC_table, gnomad_constraints, by1 = "ensembl_gene_id")


# j1 %>% filter(!is.na(gnomad_gene_length)) %>% .$hgnc_id %>% unique() %>% length()
# j2 %>% filter(!is.na(gnomad_gene_length)) %>% .$hgnc_id %>% unique() %>% length()
# HGNC_table %>% filter(!is.na(gnomad_gene_length)) %>% .$hgnc_id %>% unique() %>% length()



# get GTEX features and join with HGNC table
cat("get GTEX features...")
source("gene_score/features/scripts/gtex.R")
gtex_nTPM <- read_csv(paste0("gene_score/features/results/rna_tissue_gtex_nTPM_agg_", creation_date, ".csv"), show_col_types = FALSE)
gtex_tau <- read_csv(paste0("gene_score/features/results/rna_tissues_gtex_nTPM_agg_tau_val_", creation_date, ".csv"), show_col_types = FALSE)
# HGNC_table <- HGNC_table %>% left_join(gtex_nTPM, by = "ensembl_gene_id") %>% left_join(gtex_tau, by = "ensembl_gene_id")
HGNC_table <- left_join_rescue(HGNC_table, gtex_nTPM, by1 = "ensembl_gene_id", by2 = "symbol")
HGNC_table <- left_join_rescue(HGNC_table, gtex_tau, by1 = "ensembl_gene_id", by2 = "symbol")


# # get KidneyNetwork features and join with HGNC table
# cat("get KidneyNetwork features...")
# source("gene_score/features/scripts/kidney_network.R")
# kidney_network_sums_z_scores <- read_csv(paste0("gene_score/features/results/Kidney_Network_sums_z_scores_", creation_date, ".csv"), show_col_types = FALSE)
# HGNC_table <- HGNC_table %>% left_join(kidney_network_sums_z_scores, by = c("ensembl_gene_id" = "ensembl_id"))


# # get Nephrogenesis Atlas features and join with HGNC table
# #TODO: problem with double gene expression values for the same gene!
# cat("get Nephrogenesis Atlas features...")
# source("gene_score/features/scripts/nephrogenesis_atlas.R")
# nephrogenesis_atlas <- read_csv(paste0("gene_score/features/results/fetal_avg_expr_nephrogenesis_atlas_", creation_date, ".csv"), show_col_types = FALSE)
# HGNC_table <- HGNC_table %>% left_join(nephrogenesis_atlas, by = c("hgnc_id_int" = "hgnc_id"))


# get number of paralogues (above xth percentile Query% and Target%) and join with HGNC table
cat("get nunber of close paralogues...")
source("gene_score/features/scripts/paralogues.R")
no_paralogues <- read_csv(paste0("gene_score/features/results/paralogues_95_85_75_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(no_paralogues, by = "ensembl_gene_id")


# get promoter CpG observed-to-expected-ratio and join with HGNC table
cat("get promoter CpG observed-to-expected-ratio...")
source("gene_score/features/scripts/promoter_CpG_o2e_ratio.R")
promoter_CpG_o2e_ratio <- read_csv(paste0("gene_score/features/results/canonical_promoter_CpG_obs_to_exp_ratio_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(promoter_CpG_o2e_ratio, by = "ensembl_gene_id")


# get average exons CpG observed-to-expected-ratio and join with HGNC table
cat("get average exons CpG observed-to-expected-ratio...")
source("gene_score/features/scripts/exon_CpG_o2e_ratio.R")
exons_CpG_o2e_ratio <- read_csv(paste0("gene_score/features/results/canonical_ts_exons_CpG_obs_to_exp_ratio_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(exons_CpG_o2e_ratio, by = "ensembl_gene_id")


# get exon and promoter conservation scores and join with HGNC table
cat("get exon and promoter conservation scores...")
source("gene_score/features/scripts/exon_and_prom_conservation.R")
avg_phasCons_ex <- read_csv(paste0("gene_score/features/results/avg_phasCons_scores_per_transcript_", creation_date, ".csv"), show_col_types = FALSE)
avg_phasCons_prom <- read_csv(paste0("gene_score/features/results/avg_phasCons_promoter_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(avg_phasCons_ex, by = "ensembl_gene_id")
HGNC_table <- HGNC_table %>% left_join(avg_phasCons_prom, by = "ensembl_gene_id")


# get annotation for which genes are associated with MGI MPO MP_0005367
cat("get MGI MPO annotation...")
source("gene_score/features/scripts/mgi_mpo.R")
mgi_mpo_kidney <- read_csv(paste0("gene_score/features/results/mgi_human_genes_associated_MP_0005367_" , creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(mgi_mpo_kidney, by = "entrez_id")



# write results
write.csv(HGNC_table, paste0("gene_score/features/results/gene_features_", creation_date, ".csv"), row.names = FALSE)








write.xlsx(data.frame(features=names(HGNC_table)),
           file = "feature_names.xlsx",
           row.names = F)


data.frame(features=names(HGNC_table))



# TODO
# add descartes percent expression, tau, and kidney genes






# NOTES - TO BE DELETED!!
# 
# 
# y <- gnomad_constraints
# x <- HGNC_table
# by1 <- "ensembl_gene_id"
# 
# 
# 
# 
# left_join_rescue_symbol(df1, df2, "A", "B")
# 
# 
# # Create a sample dataframe
# df1 <- data.frame(A = c(1,2,3,4,5, 100), B = c("a", "b", "c", "d", "e", "f"), C = 5:10)
# df2 <- data.frame(A = c(1,2,3,10,6, 100), B = c("z", "b", "x", "d", "r", "l"), D = c(33, 31, 23, 45, 44, 100))
# df1
# df2
# 
# df1 %>% left_join(df2, by=c("A", "B"))
# 
# 
# left_join_rescue <- function(x, y, by_1_l, by_1_r, by_2_l, by_2_r) {
#   y_sub <- y %>% dplyr::select(-{{ by_2_r }})
#   colnames(y_sub)[which(colnames(y_sub) %in% by_1_r)] <- by_1_l
#   
#   joined_sub <- left_join(x, y_sub, by = by_1_l)
#   not_joined_ids <- setdiff(data.frame(y)[, by_1_r], data.frame(x)[, by_1_l])
# 
#   # 
#   return(not_joined_ids)
# 
# }
# 
# y <- gnomad_constraints
# x <- HGNC_table
# by_1_l <- "ensembl_gene_id"
# by_1_r <- "gene_id"
# by_2_l <- "symbol"
# by_2_r <- "gene"
# 
# j <- left_join_rescue(HGNC_table, gnomad_constraints, "ensembl_gene_id", "gene_id", "symbol", "gene")
# 
# setdiff(as.vector(HGNC_table[, by_1_l]), as.vector(gnomad_constraints[, by_1_r]))
# 
# 
# setdiff(df[, "A"], df[, "C"])
# 
# 
# setdiff(gnomad_constraints[, by_1_r], data.frame(HGNC_table)[, by_1_l])
# 
# 
# 
# 
# my_function <- function(A, dataframe) {
#   selected_column <- dataframe %>%
#     dplyr::select({{ A }})
#   
#   return(selected_column)
# }
# 
# 
# # Create a sample dataframe
# df <- data.frame(A = 1:5, B = letters[1:5], C = 5:9)
# 
# # Call the function to select the "B" column from the dataframe
# result <- my_function("B", df)
# 
# # Print the result
# print(result)
# 
# 
# library(dplyr)
# library(rlang)
# 
# my_function <- function(var1, var2, x, y_sub) {
#   joined_sub <- x %>%
#     left_join(y_sub, by = setNames(!!!syms(c(var1, var2))))
#   
#   return(joined_sub)
# }
# 
# 
# # Create sample dataframes
# x <- data.frame(var1 = c(1, 2, 3), var2 = c("A", "B", "C"), value = c(10, 20, 30))
# y_sub <- data.frame(var2 = c("A", "C"), value2 = c(100, 300))
# 
# # Call the function to perform the left join
# result <- my_function("var1", "var2", x, y_sub)
# 
# # Print the result
# print(result)
# 
# 
# 
# # SAME PROBLEM FOR GTEX => check!
# 

