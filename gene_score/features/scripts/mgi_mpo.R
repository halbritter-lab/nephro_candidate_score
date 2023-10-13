# Title: Genes associated with mouse phenotype ontology "MP:0005367" (= "renal/urinary system phenotype") and its children

# load libraries
library(tidyverse)	
library(jsonlite) 
library(config)

# read configs
config_vars <- config::get(file = "config.yml")
script_path <- "gene_score/features"

# save current working directory
wd_bef_script_exe <- getwd()

# set working directory
setwd(file.path(config_vars$PROJECT_DIR, script_path))


# download mouse phenotype ontology (MPO)
mgi_mpo_url <- "https://www.informatics.jax.org/downloads/reports/mp.json"
download.file(mgi_mpo_url, 
              destfile = paste0("raw/mgi_mpo_", config_vars$creation_date, ".json"))

# load MPO data
mp <- fromJSON(paste0("raw/mgi_mpo_", config_vars$creation_date, ".json"))

# get parent-child relations
relations <- mp$graphs$edges[[1]] %>% 
  rowwise() %>% 
  mutate(parent=str_remove(obj, "http://purl.obolibrary.org/obo/"), child=str_remove(sub, "http://purl.obolibrary.org/obo/")) %>% 
  dplyr::select(-obj, -sub, -pred)
  
# define functions to get all children of a given MPO-term 
mpo_get_direct_children <- function(parent_term, relation_df){
  direct_children <- relations %>% 
    filter(parent == parent_term) %>% 
    .$child
  return(direct_children)
}

mpo_all_children_from_term <- function(term_input, relation_df) {
  children_list <- mpo_get_direct_children(term_input, relation_df)
  mpo_all_children <<- append(mpo_all_children, term_input)
  if(length(children_list)!=0)
  {
    for (p in children_list) {
      mpo_all_children <<- append(mpo_all_children, p)
      Recall(p)
    }
  }
  return(unique(mpo_all_children))
}

# get all children of term "MP:0005367" (= "renal/urinary system phenotype")
mpo_all_children <- c()
kidney_mpo_terms <- mpo_all_children_from_term("MP_0005367", relations) %>% 
  str_replace("_", ":")


# download and load mouse/human orthology with phenotype annotations
download_url <- "https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"
download.file(url = download_url, 
              destfile = paste0("raw/HMD_HumanPhenotype_", config_vars$creation_date, ".rpt"))
hmd_hp <- read.table(paste0("raw/HMD_HumanPhenotype_", config_vars$creation_date, ".rpt"), sep="\t", header=FALSE) %>% dplyr::select(-V6)
colnames(hmd_hp) <- c("human_marker_symbol", "human_entrez_id", "mouse_marker_symbol", "mgi_marker_accession_id", "mpo_id")


# download and load genotypes and mammalian phenotype annotations
download_url <- "https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt"
download.file(url = download_url, 
              destfile = paste0("raw/MGI_PhenoGenoMP_", config_vars$creation_date, ".rpt"))

mgi_pg_mp <- read.delim(paste0("raw/MGI_PhenoGenoMP_", config_vars$creation_date,  ".rpt"), header=FALSE)
colnames(mgi_pg_mp) <- c("allellic_composition", "allele_symbols", "genetic_background", "mpo_id", "pubmed_id", "mgi_marker_accession_id")


# check which genotypes are associated with "MP:0005367" (= "renal/urinary system phenotype") or any of its children
mgi_pg_mp <- mgi_pg_mp %>% 
  rowwise() %>% 
  mutate(mpo_kid_pos = mpo_id %in% kidney_mpo_terms) 

# kidney positive genotypes 
kid_pos_gt <- mgi_pg_mp %>% 
  filter(mpo_kid_pos == TRUE)

# filter out polygenic genotypes, check for heterozygosity/homozygosity in monogenic phenotypes
kid_pos_gt <- kid_pos_gt %>% 
  rowwise() %>% 
  filter(length(strsplit(mgi_marker_accession_id, split = "|", fixed = TRUE)[[1]]) == 1) %>% # remove polygenic genotypes
  mutate(wt_alleles = sum(grepl("<\\+>", (strsplit(allele_symbols, "\\|")[[1]]))), 
         mpo_kidney_gene_het = (wt_alleles == 1),
         mpo_kidney_gene_hom = (wt_alleles == 0))


# heterozygous and homozygous genotypes associated with "MP:0005367" 
het_kid_genes <- kid_pos_gt %>% filter(mpo_kidney_gene_het == TRUE) %>% .$mgi_marker_accession_id %>% unique()
hom_kid_genes <- kid_pos_gt %>% filter(mpo_kidney_gene_hom == TRUE) %>% .$mgi_marker_accession_id %>% unique()

# annotate raw dataframe
mgi_pg_mp <- mgi_pg_mp %>% 
  dplyr::select(mgi_marker_accession_id) %>% 
  distinct() %>% 
  rowwise() %>% 
  filter(length(strsplit(mgi_marker_accession_id, split = "|", fixed = TRUE)[[1]]) == 1) %>% # remove polygenic genotypes
  mutate(mgi_kid = case_when(mgi_marker_accession_id %in% het_kid_genes ~ 2,
                             mgi_marker_accession_id %in% hom_kid_genes ~ 1,
                             .default = 0)) %>% 
  left_join(hmd_hp[, c("human_entrez_id", "human_marker_symbol", "mgi_marker_accession_id")], by = "mgi_marker_accession_id", relationship = "many-to-many") %>% # Note: some human entrez ids have multiple mouse ids
  filter(!is.na(human_entrez_id)) %>% 
  dplyr::select(entrez_id = human_entrez_id, symbol = human_marker_symbol, mgi_kid) %>% 
  distinct() %>% 
  group_by(entrez_id, symbol) %>% 
  summarize(max_mgi_kid = max(mgi_kid, na.rm = TRUE)) # if there are two values for one human entrez_id, take the highest (2 > 1 > 0, heterozygous > homozygous > no kidney phenotype)


# write results
write.csv(mgi_pg_mp, 
          paste0("results/mgi_human_genes_associated_MP_0005367_" , config_vars$creation_date, ".csv"), 
          row.names = FALSE)

# set back former working directory
setwd(wd_bef_script_exe)
