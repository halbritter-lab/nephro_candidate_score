# Title: gnomad loss of function values

# load libraries
library(tidyverse)
library(R.utils)

# download and unzip gnomad constraint data
download_url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

download.file(url = download_url, 
              destfile = "gene_score/features/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")

gunzip(filename = "gene_score/features/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", 
       destname = "gene_score/features/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt")

# extract Z- and pLI-scores
gnomad_lof <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene.txt") %>%
  dplyr::select(gene, gene_id, gene_type, syn_z, mis_z, pLI)

# write results
write.csv(gnomad_lof, paste0("gene_score/features/results/gnomad_Z_pLI_", creation_date, ".csv"), row.names = FALSE)
