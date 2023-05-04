# Title: GTEx expression values of different tissues

# load libraries
library(tidyverse)

# download and unzip GTEx RNA expression data of different tissues
gtex_download_url <- "https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip"
download.file(gtex_download_url, 
              destfile = paste0("gene_score/features/raw/rna_tissue_gtex_", creation_date, ".tsv.zip"))

unzip(zipfile = paste0("gene_score/features/raw/rna_tissue_gtex_", creation_date, ".tsv.zip"))  

# load data
rna_tissue_gtex_nTPM <- read.delim(paste0("gene_score/features/raw/rna_tissue_gtex_", creation_date, ".tsv")) %>% 
  dplyr::select(Gene, Tissue, nTPM) %>% 
  spread(key = Tissue, value = nTPM) 

# dataset contains expression values for all tissues for 19764 genes - except for retina: here for 20090 genes
# => remove genes, that only have expression values for retina, but no other tissues
rna_tissue_gtex_nTPM <- rna_tissue_gtex_nTPM[complete.cases(rna_tissue_gtex_nTPM), ]

# aggregate brain tissues
brain_regions <- c("amygdala", "caudate", "cerebellum", "cerebral cortex", "hippocampus", "hypothalamus", "nucleus accumbens", "putamen", "substantia nigra")
brain_nTPM <- rna_tissue_gtex_nTPM %>%
  dplyr::select(Gene, all_of(brain_regions))
brain_nTPM_med <- data.frame(Gene = brain_nTPM$Gene, brain_median = apply(brain_nTPM[, -1], 1, median))

# join median brain nTPM values with other tissue df
rna_tissue_gtex_nTPM_agg <- rna_tissue_gtex_nTPM %>% 
  dplyr::select(-brain_regions) %>% 
  left_join(brain_nTPM_med, by = "Gene") %>%
  column_to_rownames(var = "Gene")

# calculate the robust tissue specifity measure tau
# calculate the normalized tau values according to Yanai et al. for all genes
N <- ncol(rna_tissue_gtex_nTPM_agg)
max_row <- apply(rna_tissue_gtex_nTPM_agg, 1, max)
norm_data <- rna_tissue_gtex_nTPM_agg / max_row
tau_values <- sapply(1:nrow(norm_data), function(i) {
  xi <- norm_data[i, ]
  sum((1 - xi) / (N - 1))
})

tau_df <- data.frame(gene = rownames(norm_data), tau = tau_values)

# rownames back to column
rna_tissue_gtex_nTPM_agg <- rna_tissue_gtex_nTPM_agg %>% 
  rownames_to_column(var = "Gene") %>% 
  rename_all(~ str_replace_all(.x, " ", "_"))

# write results - nTPM values
write.csv(rna_tissue_gtex_nTPM_agg, paste0("gene_score/features/results/rna_tissue_gtex_nTPM_agg_" , creation_date, ".csv"), row.names=FALSE)

# write results - tau values
write.csv(tau_df, paste0("gene_score/features/results/rna_tissues_gtex_nTPM_agg_tau_val_" , creation_date, ".csv"), row.names=FALSE)



