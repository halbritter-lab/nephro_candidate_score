# Title: Dispensible genes - genes with homozygous loss of function variants, not associated with any phenotype in OMIM and not listed in kidney-genetics

# load libraries
library(tidyverse)
library(readr)
library(jsonlite) # required in hgnc-functions.R


# read configs
config_vars <- config::get(file = "config.yml")
script_path <- "gene_score/labels"

# save current working directory
wd_bef_script_exe <- getwd()

# set working directory
setwd(file.path(config_vars$PROJECT_DIR, script_path))

# source additional functions
source(file.path(config_vars$PROJECT_DIR, "gene_score", "hgnc_functions.R"))


# download and unzip file with homozygus loss of function variants from gnomad publication
download_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/bin/41586_2020_2308_MOESM4_ESM.zip"
download.file(url = download_url,
              destfile = paste0("raw/41586_2020_2308_MOESM4_ESM_", config_vars$creation_date, ".zip"))
unzip(zipfile = paste0("raw/41586_2020_2308_MOESM4_ESM_", config_vars$creation_date, ".zip"),
      exdir = "raw/")

# load homozygous knockout genes
hom_ko_genes <- read.table("raw/supplement/supplementary_dataset_7_hom_ko_genes.txt")

# get HGNC ID from gene symbol
hom_ko_genes  <- hgnc_id_from_symbol_grouped(tibble(value = hom_ko_genes$V1)) %>%
  tibble() %>%
  dplyr::rename(c("hgnc_id_int" = ".")) %>%
  drop_na(hgnc_id_int)




# download OMIM genemap
# OMIM links to genemap2 file needs to be set in config and applied for at https://www.omim.org/downloads
download.file(url = config_vars$omim_download_url,
              destfile = paste0("raw/genemap2_", config_vars$creation_date, ".txt"))

# load file and clean column names
names_col <- read_tsv(paste0("raw/genemap2_", config_vars$creation_date, ".txt"),
                     col_names = FALSE,
                     skip = 3,
                     n_max = 1,
                     col_types = cols())
clean_names <- gsub(" ", "_", gsub("# ", "", names_col[1, ]))

# read omim genemap
omim_genes_hg38 <-  read.delim2(paste0("raw/genemap2_", config_vars$creation_date, ".txt"), 
                                header = FALSE, 
                                comment.char = "#") %>% 
  dplyr::mutate_all(~ ifelse(. == "", NA, .))
colnames(omim_genes_hg38) <- clean_names

# filter for morbid genes (genes associated with a phenotype in OMIM)
# omim_hg38_morbid <- omim_genes_hg38 %>% 
#   filter(!is.na(Phenotypes), !is.na(Approved_Gene_Symbol)) # Note 2023-10-12: there are no genes in ho_ko_genes that have a Phenotypes but no Approved_Symbol

# filter for morbid genes (genes associated with a phenotype in OMIM)
omim_hg38_morbid <- omim_genes_hg38 %>% 
  filter(!is.na(Phenotypes))

# morbid genes with 'Approved_Gene_Symbol'
omim_hg38_morbid_with_approved_symbol <- omim_hg38_morbid %>% filter(!is.na(Approved_Gene_Symbol))

# get hgnc ids of morbid genes with 'Approved_Gene_Symbol'
omim_hg38_morbid_was_hgnc <- hgnc_id_from_symbol_grouped(tibble(value = omim_hg38_morbid_with_approved_symbol$Approved_Gene_Symbol)) %>%
  tibble() %>% 
  dplyr::rename(c("hgnc_id_int" = ".")) %>% 
  drop_na(hgnc_id_int)

# morbid genes without 'Approved_Gene_Symbol'
omim_hg38_morbid_without_approved_symbol <- omim_hg38_morbid %>% filter(is.na(Approved_Gene_Symbol))

# get hgnc ids of morbid genes without 'Approved_Gene_Symbol' (long processing time)
omim_hg38_morbid_without_approved_symbol <- separate_rows(omim_hg38_morbid_without_approved_symbol, `Gene/Locus_And_Other_Related_Symbols`, sep = ", ")
omim_hg38_morbid_woas_hgnc <- hgnc_id_from_symbol_grouped(tibble(value = unique(omim_hg38_morbid_without_approved_symbol$`Gene/Locus_And_Other_Related_Symbols`))) %>%
  tibble() %>% 
  dplyr::rename(c("hgnc_id_int" = ".")) %>% 
  drop_na(hgnc_id_int)

# combine morbid genes with and without 'Approved_Gene_Symbol'
omim_hg38_morbid_hgnc <- c(omim_hg38_morbid_was_hgnc$hgnc_id_int, omim_hg38_morbid_woas_hgnc$hgnc_id_int)

# get dispensible genes (homozygous ko genes without 'Phenotypes' entry in OMIM and not listed in kidney-genetics)
disp_genes <- hom_ko_genes %>% 
  filter(!(hgnc_id_int %in% unique(omim_hg38_morbid_hgnc)), !(hgnc_id_int %in% kid_gen$hgnc_id_int)) %>% 
  distinct()

# write results
write.csv(disp_genes, paste0("results/dispensible_genes_", config_vars$creation_date, ".csv"), row.names = FALSE)

# set back former working directory
setwd(wd_bef_script_exe)
