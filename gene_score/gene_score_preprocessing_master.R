# Title: Gene score data preprocessing master script

# load libraries
library(tidyverse)
library(utils)
library(R.utils)
library(config)
library(jsonlite)

# define relative script path
project_topic <- "nephrology"
project_name <- "nephro_candidate_score"
script_path <- "/gene_score/"

# read configs
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"),
                           config = project_topic)

# set working directory
setwd(paste0(config_vars$projectsdir, project_name, script_path))

# store wd
init_wd <- getwd()

# source additional functions
source("hgnc_functions.R")
source("helper_functions.R")

##### HGNC table ##### 
# download HGNC gene table from github repository "kidney-genetics"
gene_table_url <- paste0("https://github.com/halbritter-lab/kidney-genetics/blob/main/analyses/B_AnnotationHGNC/results/non_alt_loci_set_coordinates.", config_vars$hgnc_gt_version, ".csv.gz?raw=TRUE")
download.file(url = gene_table_url, destfile = paste0("raw/HGNC_", config_vars$hgnc_gt_version, ".csv.gz"))

# load HGNC gene table and filter for protein-coding genes # TODO: change paths
HGNC_table <- read_csv(paste0("raw/HGNC_", config_vars$hgnc_gt_version, ".csv.gz"), 
                       show_col_types = FALSE) %>% 
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
setwd(init_wd)

# get dispensible genes #TODO: test with stable internet connection
cat("get dispensible genes...")
source("labels/scripts/dispensable_genes.R")
setwd(init_wd)


##### FEATURES ##### 
# get cellXgene features and join with HGNC table
cat("get cellXgene features...")
source("features/scripts/cellxgene.R")
setwd(init_wd)
cellXgene_ds1 <- read_csv(paste0("features/results/cellxgene_expr_0b4a15a7-4e9e-4555-9733-2423e5c66469_", config_vars$creation_date, ".csv.gz"), 
                          show_col_types = FALSE)
HGNC_table <- left_join_rescue_symbol(HGNC_table, cellXgene_ds1, by1 = "ensembl_gene_id") 


# get Descartes features and join with HGNC table
cat("get descartes fetal kidney features...")
source("features/scripts/descartes.R")
setwd(init_wd)

descartes_fetal_kid_tau <- read_csv(paste0("features/results/descartes_fetal_kidney_tau_", config_vars$creation_date, ".csv.gz"), 
                                    show_col_types = FALSE)
descartes_fetal_kid_pe <- read_csv(paste0("features/results/descartes_fetal_kidney_percent_expression_", config_vars$creation_date, ".csv.gz"), 
                                   show_col_types = FALSE)
descartes_fetal_kid_nptm <- read_csv(paste0("features/results/descartes_fetal_nptm_", config_vars$creation_date, ".csv.gz"), 
                                     show_col_types = FALSE)

HGNC_table <- left_join_rescue_symbol(HGNC_table, descartes_fetal_kid_tau, by1 = "ensembl_gene_id") 
HGNC_table <- left_join_rescue_symbol(HGNC_table, descartes_fetal_kid_pe, by1 = "ensembl_gene_id")
HGNC_table <- left_join_rescue_symbol(HGNC_table, descartes_fetal_kid_nptm, by1 = "ensembl_gene_id")
  

# get gnomAD features and join with HGNC table
cat("get gnomad features...")
source("features/scripts/gnomad.R")
setwd(init_wd)

gnomad_constraints <- read_csv(paste0("features/results/gnomad_constraints_", config_vars$creation_date, ".csv.gz"), 
                               show_col_types = FALSE)
HGNC_table <- left_join_rescue_symbol(HGNC_table, gnomad_constraints, by1 = "ensembl_gene_id")


# get GTEX features and join with HGNC table
cat("get GTEX features...")
source("features/scripts/gtex.R")
setwd(init_wd)

gtex_nTPM <- read_csv(paste0("features/results/rna_tissue_gtex_nTPM_agg_", config_vars$creation_date, ".csv.gz"), 
                      show_col_types = FALSE)
gtex_tau <- read_csv(paste0("features/results/rna_tissues_gtex_nTPM_agg_tau_val_", config_vars$creation_date, ".csv.gz"), 
                     show_col_types = FALSE)
HGNC_table <- left_join_rescue_symbol(HGNC_table, gtex_nTPM, by1 = "ensembl_gene_id")
HGNC_table <- left_join_rescue_symbol(HGNC_table, gtex_tau, by1 = "ensembl_gene_id")


# get number of close paralogues (above xth percentile Query% and Target%) and join with HGNC table
cat("get nunber of close paralogues...")
source("features/scripts/paralogues.R")
setwd(init_wd)

no_paralogues <- read_csv(paste0("features/results/paralogues_95_85_75_", config_vars$creation_date, ".csv.gz"), 
                          show_col_types = FALSE)
HGNC_table <- left_join_rescue_symbol(HGNC_table, no_paralogues, by1 = "ensembl_gene_id")


# get promoter CpG observed-to-expected-ratio and join with HGNC table
cat("get promoter CpG observed-to-expected-ratio...")
source("features/scripts/promoter_CpG_o2e_ratio.R")
setwd(init_wd)

promoter_CpG_o2e_ratio <- read_csv(paste0("features/results/canonical_promoter_CpG_obs_to_exp_ratio_", config_vars$creation_date, ".csv.gz"), 
                                   show_col_types = FALSE)
HGNC_table <- left_join_rescue_symbol(HGNC_table, promoter_CpG_o2e_ratio, by1 = "ensembl_gene_id")


# get average exons CpG observed-to-expected-ratio and join with HGNC table
cat("get average exons CpG observed-to-expected-ratio...")
source("features/scripts/exon_CpG_o2e_ratio.R")
setwd(init_wd)

exons_CpG_o2e_ratio <- read_csv(paste0("features/results/canonical_ts_exons_CpG_obs_to_exp_ratio_", config_vars$creation_date, ".csv.gz"), 
                                show_col_types = FALSE)
HGNC_table <- left_join_rescue_symbol(HGNC_table, exons_CpG_o2e_ratio, by1 = "ensembl_gene_id")


# get exon and promoter conservation scores and join with HGNC table
cat("get exon and promoter conservation scores...")
source("features/scripts/exon_and_prom_conservation.R")
setwd(init_wd)

avg_phasCons_ex <- read_csv(paste0("features/results/avg_phasCons_scores_per_transcript_", config_vars$creation_date, ".csv.gz"), 
                            show_col_types = FALSE)
avg_phasCons_prom <- read_csv(paste0("features/results/avg_phasCons_promoter_", config_vars$creation_date, ".csv.gz"), 
                              show_col_types = FALSE)
HGNC_table <- left_join_rescue_symbol(HGNC_table, avg_phasCons_ex, by1 = "ensembl_gene_id")
HGNC_table <- left_join_rescue_symbol(HGNC_table, avg_phasCons_prom, by1 = "ensembl_gene_id")


# get annotation for which genes are associated with MGI MPO MP_0005367
cat("get MGI MPO annotation...")
source("features/scripts/mgi_mpo.R")
setwd(init_wd)

mgi_mpo_kidney <- read_csv(paste0("features/results/mgi_human_genes_associated_MP_0005367_" , config_vars$creation_date, ".csv.gz"), 
                           show_col_types = FALSE)
HGNC_table <- left_join_rescue_symbol(HGNC_table, mgi_mpo_kidney, by1 = "entrez_id")


# write results
all_gene_features <- HGNC_table %>% 
  dplyr::select(-hgnc_id, -entrez_id, -ensembl_gene_id, -symbol, -alias_symbol, -prev_symbol)

write.csv(all_gene_features, 
          paste0("features/results/gene_features_", config_vars$creation_date, ".csv"), 
          row.names = FALSE)

gzip(paste0("features/results/gene_features_", config_vars$creation_date, ".csv"),
     overwrite = TRUE)


