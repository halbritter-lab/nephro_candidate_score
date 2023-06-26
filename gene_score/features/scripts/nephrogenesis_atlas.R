# Title: Average ssRNA expression values in kidney genesis from "Human Nephrogenesis Atlas" (https://sckidney.flatironinstitute.org)

# load libraries
library(dplyr)
library(Seurat)
library(data.table)

# download Seurat object
nephrogenesis_atlas_url <- "https://sckidney.flatironinstitute.org/media/Human_nephrogenesis_atlas.Robj"
download.file(nephrogenesis_atlas_url,
              destfile = paste0("gene_score/features/raw/Human_nephrogenesis_atlas", creation_date, ".Robj"))

# load data
load(paste0("gene_score/features/raw/Human_nephrogenesis_atlas_", creation_date, ".Robj"))
fetal <- Sharing_is_Caring

# map cluster values to cluster names
fetal@active.ident <- plyr::mapvalues(
  x = fetal@active.ident,
  from = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"),
  to = c("Nephron_progen_1", 
         "Nephron_progen_2", 
         "Committed_progen_1", 
         "Committed_progen_2", 
         "Pretub_aggreg_1", 
         "Pretub_aggreg_2", 
         "Prox_RV", 
         "Podo_of_SSB", 
         "Parietal_epithel_of_SSB", 
         "Cycling_prox_early_nephron", 
         "Prox_medial_RV", 
         "PT_SSB", 
         "Distal_pretub_aggregate_RV", 
         "Distal_RV_CSB", 
         "Cycling_distal_early_nephron",
         "Dist_tub_SSB", 
         "CNT_SSB", 
         "LOH_MD_SSB")
  )
Idents(fetal) <- fetal@active.ident

# normalize data
fetal <- NormalizeData(fetal, 
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)

# get cluster average expression
cluster_avg <- AverageExpression(fetal, return.seurat = FALSE, 
                                      slot = "data", 
                                      verbose = TRUE) %>% .$RNA %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene")

# prefilter for protein coding genes from HGNC table
all_prot_coding_gene_symbols_cap <- toupper(all_prot_coding_gene_symbols)

cluster_avg <- cluster_avg %>% 
  rowwise %>% 
  mutate(gene = toupper(gene)) %>% 
  filter(gene %in% all_prot_coding_gene_symbols_cap)

# annotate with HGNC IDs
cluster_avg$hgnc_id <- hgnc_id_from_symbol_grouped(tibble(value = cluster_avg$gene)) 

cluster_avg <- cluster_avg %>% dplyr::select(-gene)


# write results
write.csv(cluster_avg, 
          paste0("gene_score/features/results/fetal_avg_expr_nephrogenesis_atlas_", creation_date, ".csv")) 
