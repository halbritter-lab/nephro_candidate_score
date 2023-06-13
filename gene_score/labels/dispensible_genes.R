# UNFINISHED - CHECK TODOs!
# General TODO: shorten
# TODO: check ALPK1 in final df

# Title: Dispensible genes - genes with homozygous loss of fucnction variants, but not associated with any phenotype in OMIM

# load libraries
library(tidyverse)
library(readr)

# download and unzip file with homozygous loss of function variants from gnomad publication
download_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/bin/41586_2020_2308_MOESM4_ESM.zip"
download.file(url = download_url, 
              destfile = paste0("gene_score/labels/raw/41586_2020_2308_MOESM4_ESM_", creation_date, ".zip"))
unzip(paste0("gene_score/labels/raw/41586_2020_2308_MOESM4_ESM_", creation_date, ".zip"))

# load data
hom_lof_genes <- read.table(paste0("gene_score/labels/raw/supplement/supplementary_dataset_7_hom_ko_genes.txt"))

# TODO: add convert function!!
hom_lof_genes_converted <- convert_genes_to_approved(hom_ko_genes_leitao$X1)
hom_lof_genes_final <- hom_ko_genes_converted[[1]]$approved_symbol



# TODO: add note in master file to download of genemap2.txt from OMIM manually (registration required, access expires after 1 year)
names_col <- read_tsv(paste0("gene_score/labels/raw/genemap2", creation_date, ".txt"),
                     col_names = FALSE,
                     skip = 3,
                     n_max = 1,
                     col_types = cols())

omim_genes_hg38 <-  read_tsv(paste0("gene_score/labels/raw/genemap2", creation_date, ".txt"),
                           skip = 2,
                           col_names = TRUE,
                           comment = "#",
                           col_types = cols()) 

colnames(omim_genes_hg38) <-  names_col

# annotate with "morbid" if gene is associated with a phenotype in OMIM
omim_genes_hg38_annot <-  omim_genes_hg38 %>%
  rename(Chromosome = "# Chromosome") %>%
  rowwise() %>%
  mutate(morbid = case_when(!is.na(`Approved Symbol`) & !is.na(Phenotypes) ~ T, T ~ F))

# OMIM morbid genes (subset of genes that have "Approved Symbol" & "Phenotype")
omim_morbid_all <– omim_genes_hg38_annot %>% 
  filter(morbid == T)

# dispensible genes
# TODO: replace HGNC_symbols_allChr_approved with our df
# TODO: filter out positive genes
disp_df <- HGNC_symbols_allChr_approved %>% 
  filter(!(`Approved symbol` %in% omim_morbid_all$`Approved Symbol`), `Approved symbol` %in% hom_ko_genes_final) 
