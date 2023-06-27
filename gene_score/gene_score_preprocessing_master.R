# Title: Gene score data preprocessing master script

library(tidyverse)
library(utils)
library(R.utils)
source("config.R")

# attach tidyr to second position for resolving the select() problem
if ("package:tidyr" %in% search()) detach("package:tidyr", unload = TRUE)
attachNamespace("tidyr", pos = 2)
library(dplyr)

# TODO: missing value imputation


# download HGNC gene table from github repository "kidney-genetics"
# hgnc_gt_version <- "2023-06-21"
gene_table_url <- paste0("https://github.com/halbritter-lab/kidney-genetics/raw/main/analyses/A_AnnotationHGNC/results/non_alt_loci_set_coordinates.", hgnc_gt_version, ".csv.gz")
download.file(url = gene_table_url, destfile = paste0("gene_score/raw/HGNC_", hgnc_gt_version, ".csv.gz"))
gunzip(filename = paste0("gene_score/raw/HGNC_", hgnc_gt_version, ".csv.gz"), 
       destname = paste0("gene_score/raw/HGNC_", hgnc_gt_version, ".csv"))

# load HGNC gene table and filter for protein-coding genes # TODO: change paths
HGNC_table <- read.csv("gene_score/raw/HGNC_2023-06-21.csv") %>% 
  filter(locus_group == "protein-coding gene") %>% 
  dplyr::select(hgnc_id, entrez_id, ensembl_gene_id, symbol, alias_symbol, prev_symbol) %>% 
  rowwise() %>% 
  mutate(hgnc_id_int = as.integer(str_replace(hgnc_id, "HGNC:", "")))

# extract all protein coding gene symbols
all_prot_coding_gene_symbols <- unlist(strsplit(c(HGNC_table$symbol, HGNC_table$alias_symbol, HGNC_table$prev_symbol), "\\|")) %>% unique()


##### LABELS ##### 
# get positive genes
cat("get positive genes...")
source("gene_score/labels/positive_genes.R")

# get dipensible genes #TODO: test with stable internet connection
cat("get dispensible genes...")
source("gene_score/labels/dispensible_genes.R")


##### FEATURES ##### 
# get cellXgene features and join with HGNC table
cat("get cellXgene features...")
source("gene_score/features/scripts/cellxgene.R")
cellXgene_ds1 <- read_csv(paste0("gene_score/features/results/cellxgene_expr_0b4a15a7-4e9e-4555-9733-2423e5c66469_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(cellXgene_ds1, by = c("ensembl_gene_id" = "ensembl_id"))


# get gnomAD features and join with HGNC table
cat("get gnomad features...")
source("gene_score/features/scripts/gnomad.R")
gnomad_constraints <- read_csv(paste0("gene_score/features/results/gnomad_constraints_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(gnomad_constraints, by = c("ensembl_gene_id" = "gene_id"))

# get GTEX features and join with HGNC table
cat("get GTEX features...")
source("gene_score/features/scripts/gtex.R")
gtex_nTPM <- read_csv(paste0("gene_score/features/results/rna_tissue_gtex_nTPM_agg_", creation_date, ".csv"), show_col_types = FALSE)
gtex_tau <- read_csv(paste0("gene_score/features/results/rna_tissues_gtex_nTPM_agg_tau_val_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(gtex_nTPM, by = c("ensembl_gene_id" = "Gene")) %>% left_join(gtex_tau, by = c("ensembl_gene_id" = "gene"))


# get KidneyNetwork features and join with HGNC table
cat("get KidneyNetwork features...")
source("gene_score/features/scripts/kidney_network.R")
kidney_network_sums_z_scores <- read_csv(paste0("gene_score/features/results/Kidney_Network_sums_z_scores_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(kidney_network_sums_z_scores, by = c("ensembl_gene_id" = "ensembl_id"))


# get Nephrogenesis Atlas features and join with HGNC table
#TODO: problem with double gene expression values for the same gene!
cat("get Nephrogenesis Atlas features...")
source("gene_score/features/scripts/nephrogenesis_atlas.R")
nephrogenesis_atlas <- read_csv(paste0("gene_score/features/results/fetal_avg_expr_nephrogenesis_atlas_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(nephrogenesis_atlas, by = c("hgnc_id_int" = "hgnc_id"))


# get number of paralogues (above 95th percentile Query% and Target%) and join with HGNC table
cat("get nunber of close paralogues...")
source("gene_score/features/scripts/paralogues.R")
paralogues_95 <- read_csv(paste0("gene_score/features/results/paralogues_95_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(paralogues_95, by = c("ensembl_gene_id" = "ensembl_gene_id"))


# get promoter CpG observed-to-expected-ratio and join with HGNC table
cat("get promoter CpG observed-to-expected-ratio...")
source("gene_score/features/scripts/promoter_CpG_o2e_ratio.R")
promoter_CpG_o2e_ratio <- read_csv(paste0("gene_score/features/results/canonical_promoter_CpG_obs_to_exp_ratio_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(promoter_CpG_o2e_ratio, by = c("ensembl_gene_id" = "ensembl_gene_id"))


# get exon and promoter conservation scores and join with HGNC table
cat("get exon and promoter conservation scores...")
source("gene_score/features/scripts/exon_and_prom_conservation.R")
avg_phasCons_ex <- read_csv(paste0("gene_score/features/results/avg_phasCons_scores_per_transcript_", creation_date, ".csv"), show_col_types = FALSE)
avg_phasCons_prom <- read_csv(paste0("gene_score/features/results/avg_phasCons_promoter_", creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(avg_phasCons_ex, by = c("ensembl_gene_id" = "ensembl_gene_id"))
HGNC_table <- HGNC_table %>% left_join(avg_phasCons_prom, by = c("ensembl_gene_id" = "ensembl_gene_id"))


# get annotation for which genes are associated with MGI MPO MP_0005367
cat("get MGI MPO annotation...")
source("gene_score/features/scripts/mgi_mpo.R")
mgi_mpo_kidney <- read_csv(paste0("gene_score/features/results/mgi_human_genes_associated_MP_0005367_" , creation_date, ".csv"), show_col_types = FALSE)
HGNC_table <- HGNC_table %>% left_join(mgi_mpo_kidney, by = c("entrez_id" = "human_entrez_id"))




