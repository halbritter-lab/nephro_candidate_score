# UNFINISHED - CHECK TODOs!
# General TODO: shorten
# TODO: check ALPK1 in final df

# Title: Dispensible genes - genes with homozygous loss of fucnction variants, but not associated with any phenotype in OMIM

# load libraries
library(tidyverse)
library(readr)
source("../hgnc-functions.R") # TODO:  change PATH and change 'select' into 'dplyr::select' in hgnc-functions.R

# download and unzip file with homozygous loss of function variants from gnomad publication
download_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/bin/41586_2020_2308_MOESM4_ESM.zip"
download.file(url = download_url, 
              destfile = paste0("gene_score/labels/raw/41586_2020_2308_MOESM4_ESM_", creation_date, ".zip"))
unzip(paste0("gene_score/labels/raw/41586_2020_2308_MOESM4_ESM_", creation_date, ".zip"))

# load homozygous knockout genes
hom_ko_genes <- read.table(paste0("gene_score/labels/raw/supplement/supplementary_dataset_7_hom_ko_genes.txt")) %>% 
  hgnc_id_from_symbol_grouped(tibble(value = hom_lof_genes$V1)) %>% 
  tibble() %>% 
  dplyr::rename(c("hgnc_id" = ".")) %>% 
  drop_na(hgnc_id)

# download OMIM genemap
download.file(url = omim_download_url,
              destfile = paste0("gene_score/labels/raw/genemap2_", creation_date, ".txt"))

names_col <- read_tsv(paste0("gene_score/labels/raw/genemap2", creation_date, ".txt"),
                     col_names = FALSE,
                     skip = 3,
                     n_max = 1,
                     col_types = cols())
clean_names <- gsub(" ", "_", gsub("# ", "", names_col[1,]))


omim_genes_hg38 <-  read.delim2(paste0("gene_score/labels/raw/genemap2", creation_date, ".txt"), 
                                header = FALSE, 
                                comment.char = "#") %>% 
  dplyr::mutate_all(~ifelse(. == "", NA, .))

colnames(omim_genes_hg38) <-  clean_names

# filter for morbid genes (a gene associated with a phenotype in OMIM)
omim_hg38_morbid <- omim_genes_hg38 %>% filter(!is.na(Phenotypes))

# get dispensible genes
# TODO: filter out positive genes
disp_genes <- hom_ko_genes %>% 
  filter(!(`Approved symbol` %in% omim_morbid_all$`Approved Symbol`), `Approved symbol` %in% hom_ko_genes_final) 
