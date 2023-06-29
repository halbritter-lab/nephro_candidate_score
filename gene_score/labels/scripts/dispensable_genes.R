# Title: Dispensible genes - genes with homozygous loss of function variants, not associated with any phenotype in OMIM and not listed in kidney-genetics

# load libraries
library(tidyverse)
library(readr)
library(jsonlite) # required in hgnc-functions.R
source("https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/functions/hgnc-functions.R") # TODO: maybe change this

# download and unzip file with homozygus loss of function variants from gnomad publication
download_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/bin/41586_2020_2308_MOESM4_ESM.zip"
download.file(url = download_url,
              destfile = paste0("gene_score/labels/raw/41586_2020_2308_MOESM4_ESM_", creation_date, ".zip"))
unzip(zipfile = paste0("gene_score/labels/raw/41586_2020_2308_MOESM4_ESM_", creation_date, ".zip"),
      exdir = "gene_score/labels/raw/")

# load homozygous knockout genes
hom_ko_genes <- read.table("gene_score/labels/raw/supplement/supplementary_dataset_7_hom_ko_genes.txt")

# get HGNC ID from gene symbol
hom_ko_genes  <- hgnc_id_from_symbol_grouped(tibble(value = hom_ko_genes$V1)) %>%
  tibble() %>%
  dplyr::rename(c("hgnc_id" = ".")) %>%
  drop_na(hgnc_id)

# download OMIM genemap
# download.file(url = omim_download_url,
#               destfile = paste0("gene_score/labels/raw/genemap2_", creation_date, ".txt"))
response <- GET(omim_download_url)
writeBin(content(response, "raw"), paste0("gene_score/labels/raw/genemap2_", creation_date, ".txt"))

# extract and clean column names
names_col <- read_tsv(paste0("gene_score/labels/raw/genemap2_", creation_date, ".txt"),
                     col_names = FALSE,
                     skip = 3,
                     n_max = 1,
                     col_types = cols())
clean_names <- gsub(" ", "_", gsub("# ", "", names_col[1, ]))

# read omim genemap
omim_genes_hg38 <-  read.delim2(paste0("gene_score/labels/raw/genemap2_", creation_date, ".txt"), 
                                header = FALSE, 
                                comment.char = "#") %>% 
  dplyr::mutate_all(~ ifelse(. == "", NA, .))
colnames(omim_genes_hg38) <- clean_names

# filter for morbid genes (genes associated with a phenotype in OMIM)
omim_hg38_morbid <- omim_genes_hg38 %>% 
  filter(!is.na(Phenotypes), !is.na(Approved_Gene_Symbol)) # so far no genes with Phenotypes but without Approved_Symbol exists in hom_ko_genes

# get HGNC ID for omim morbid genes
omim_hg38_morbid_hgnc <- hgnc_id_from_symbol_grouped(tibble(value = omim_hg38_morbid$Approved_Gene_Symbol)) %>%
  tibble() %>% 
  dplyr::rename(c("hgnc_id" = ".")) %>% 
  drop_na(hgnc_id)

# get dispensible genes (homozygous ko genes without Phenotypes entry in OMIM and not listed in kidney-genetics)
disp_genes <- hom_ko_genes %>% 
  filter(!(hgnc_id %in% omim_hg38_morbid_hgnc$hgnc_id), !(hgnc_id %in% kid_gen$hgnc_id)) 

