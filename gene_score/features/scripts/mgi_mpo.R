# load libraries
library(tidyverse)	
library(jsonlite) 

# download mouse phenotype ontology
mgi_mpo_url <- "https://www.informatics.jax.org/downloads/reports/mp.json"
download.file(mgi_mpo_url, 
              destfile = paste0("gene_score/features/raw/mgi_mpo_", creation_date, ".json"))

# load MPO data
mp <- fromJSON(paste0("gene_score/features/raw/mgi_mpo_", creation_date, ".json"))

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

# get all children of "MP:0005367" (= "renal/urinary system phenotype")
mpo_all_children <- c()
kidney_mpo_terms <- mpo_all_children_from_term("MP_0005367", relations) %>% 
  str_replace("_", ":")

######## load mouse/human orthology with phenotype annotations
download_url <- "https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"
# download.file(url = download_url, destfile = paste0("./HMD_HumanPhenotype_", creation_date, ".rpt"))
hmd_hp <- read.table(paste0("./HMD_HumanPhenotype_", creation_date, ".rpt"), sep="\t", header=FALSE) %>% select(-V6)
#hmd_hp <- read.table(paste0("./HMD_HumanPhenotype_2023-03-02.rpt"), sep="\t", header=FALSE) %>% select(-V6)
colnames(hmd_hp) <- c("human_marker_symbol", "human_entrez_id", "mouse_marker_symbol", "mgi_marker_accession_id", "mpo_id")

######## load genotypes and mammalian phenotype annotations
download_url <- "https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt"
#download.file(url = download_url, destfile = paste0("./MGI_PhenoGenoMP_", creation_date, ".rpt"))
mgi_pg_mp <- read.delim(paste0("./MGI_PhenoGenoMP_", creation_date,  ".rpt"), header=FALSE)
# mgi_pg_mp <- read.delim("./MGI_PhenoGenoMP_2023-03-02.rpt", header=FALSE)
colnames(mgi_pg_mp) <- c("allellic_composition", "allele_symbols", "genetic_background", "mpo_id", "pubmed_id", "mgi_marker_accession_id")

####### filter genotypes for association with "MP:0005367" "renal/urinary system phenotype" and its children
mgi_pg_mp <- mgi_pg_mp %>% 
  filter(mpo_id %in% kidney_mpo_terms) 

# separate multigenic genotypes
mgi_pg_mp <- separate(data = mgi_pg_mp, col = mgi_marker_accession_id, into = c("mgi_marker_accession_id1", "mgi_marker_accession_id2"), sep = "\\|")

# get human entrez IDs correlates that are associated with "MP:0005367" and its children in mice
hmd_hp_mpo_kidney <- hmd_hp %>% 
  filter(mgi_marker_accession_id %in% unique(c(mgi_pg_mp$mgi_marker_accession_id1, mgi_pg_mp$mgi_marker_accession_id2))) %>% 
  select(human_entrez_id) %>% 
  unique %>% 
  mutate(mpo_kidney_gene = TRUE) 

  

# write results
write.csv(hmd_hp_mpo_kidney, paste0("./mgi_human_genes_associated_MP_0005367_" , creation_date, ".csv"), row.names=FALSE)


################################################################
# TODO: discuss with Bernt whether only phenotypes directly attributed to mutations/alleles should be used or also multigenic phenotypes
# Example: MGI:1339949 Adamts4
hmd_hp_mpo_kidney %>% filter(mgi_marker_accession_id == "MGI:1339949") 
# => does not include a kidney mpo term, only attributed to kidney phenotype in multigenic phenotype Adamts1<tm1Mapr>|Adamts4<tm1.1Boer>








