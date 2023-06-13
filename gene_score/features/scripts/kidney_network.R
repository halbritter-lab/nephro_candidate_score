# Title: Summed z-scores from KidneyNetwork (https://kidney.genenetwork.nl)

# load libraries and source helper function
library(tidyverse)
library(R.utils)

options(timeout=300) #download timeout = 5min

source("gene_score/features/helper_functions.R") 

# download and unzip gene-pathway prediction z-scores of bonferroni corrected significantly predicted pathways (approx. 1GB)
download_url <- "https://molgenis26.gcc.rug.nl/downloads/KidneyNetwork/gene_pathway_scores_bonf.txt.gz"
download.file(url = download_url, destfile = paste0("gene_score/features/raw/gene_pathway_scores_bonf_", creation_date, ".txt.gz"))
gunzip(filename = paste0("gene_score/features/raw/gene_pathway_scores_bonf_", creation_date, ".txt.gz"), destname = paste0("gene_score/features/raw/gene_pathway_scores_bonf_", creation_date, ".txt"))

# load KidneyNetwork prediction z-scores
kidnet <- read_delim(paste0("gene_score/features/raw/gene_pathway_scores_bonf_", creation_date, ".txt"), "\t") %>% rename(ensembl_id = "-") 

# function for summing up positive (negative z-scores) of all children of each of the following HPO-terms
sum_KidNet_z_scores <- function(hpo_terms){
  term <- str_remove(deparse(substitute(hpo_terms)), "all_children_")
  kidnet_sub <- kidnet[, c("ensembl_id", names(kidnet)[names(kidnet) %in% hpo_terms])]
  sub_joined <- data.frame(ensembl_id = kidnet_sub$ensembl_id, 
                           pos_z_sum = apply(kidnet_sub[,2:ncol(kidnet_sub)], 1, function(x) sum(x[x > 0], na.rm = TRUE)),
                           neg_z_sum = apply(kidnet_sub[,2:ncol(kidnet_sub)], 1, function(x) sum(x[x < 0], na.rm = TRUE)))
  colnames(sub_joined)[2:3] <- paste0(term, "_", colnames(sub_joined)[2:3])
  return(sub_joined)
}

# sum up positive and negative z-scores of all children of each of the following HPO-terms:
# 1. "Abnormal renal morphology", HP:0012210
all_children_list <- list()
HPO_all_children_from_term("HP:0012210")
all_children_HP0012210 <- all_children_list %>% unlist() %>% unique()
sums_z_scores_HP0012210 <- sum_KidNet_z_scores(all_children_HP0012210)

# 2. "Abnormal renal physiology", HP:0012211
all_children_list <- list()
HPO_all_children_from_term("HP:0012211")
all_children_HP0012211 <- all_children_list %>% unlist() %>% unique()
sums_z_scores_HP0012211 <- sum_KidNet_z_scores(all_children_HP0012211)

# 3. "Abnormality of the urinary system", HP:0000079 (includes the upper two HPO-terms)
all_children_list <- list()
HPO_all_children_from_term("HP:0000079")
all_children_HP0000079 <- all_children_list %>% unlist() %>% unique()
sums_z_scores_HP0000079 <- sum_KidNet_z_scores(all_children_HP0000079)

# join dataframes
sums_z_scores <- sums_z_scores_HP0000079 %>% 
  left_join(sums_z_scores_HP0012211, by="ensembl_id") %>% 
  left_join(sums_z_scores_HP0012210, by="ensembl_id") 
  
# write results
write.csv(sums_z_scores, paste0("Kidney_Network_sums_z_scores_" , creation_date, ".csv"), row.names=FALSE)


                                             
                                             
                                             
                                             
                                             

gtex_download_url <- "https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip"
download.file(gtex_download_url, 
              destfile = paste0("gene_score/features/raw/rna_tissue_gtex_", creation_date, ".tsv.zip"))
