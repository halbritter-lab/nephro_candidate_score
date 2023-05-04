# Title: Paralogues from Ensembl

# load libraries 
library(tidyverse)
library(biomaRt)

# "Paralogues are defined in Ensembl as genes for which the most common ancestor node is a duplication event
# These ancestral duplications are represented by red nodes in the gene trees.
# The table shows the taxonomic level of the ancestor duplication node,
# the Ensembl gene ID and name, the location of the paralogue, 
# and the percent of identical amino acids in the paralogue compared with the gene of interest (Target %ID)
# The identity of the gene of interest when compared with the paralogue is the query %ID.", Ensembl website

# download paralogues from Ensembl 
bm_version <- 109  # Ensembl Biomart Genes version

ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      version = bm_version)

paralogues <- getBM(attributes = c("ensembl_gene_id",
                                   "chromosome_name",
                                   "hsapiens_paralog_associated_gene_name",
                                   "hsapiens_paralog_chromosome",
                                   "hsapiens_paralog_orthology_type",
                                   "hsapiens_paralog_subtype",
                                   "hsapiens_paralog_perc_id",
                                   "hsapiens_paralog_perc_id_r1"),
                    filter = "transcript_biotype", 
                    values = "protein_coding",
                    mart = ensembl) %>%
  filter(chromosome_name %in% c(as.character(seq(1:22)), "X", "Y")) %>% 
  dplyr::select(-chromosome_name)

# Calculate 95th percentiles of the Target% and Query%
target_95_perc <-  quantile(paralogues$hsapiens_paralog_perc_id, probs = seq(0, 1, 0.01), na.rm = TRUE)[96] # 95% percentile Target%
query_95_perc <-  quantile(paralogues$hsapiens_paralog_perc_id_r1, probs = seq(0, 1, 0.01), na.rm = TRUE)[96] # 95% percentile Query%

# determine number of paralogues above the 95th percentile (target% and query%) per gene 
paralogues_95 <- paralogues %>%
  filter(hsapiens_paralog_perc_id > target_95_perc & hsapiens_paralog_perc_id_r1 > query_95_perc) %>% 
  unique() %>%
  group_by(ensembl_gene_id) %>%
  summarize(no_paralogues_95 = n())

# write results
write.csv(paralogues_95, 
          paste0("gene_score/features_results/paralogues_95_", creation_date, ".csv"), 
          row.names = FALSE)


