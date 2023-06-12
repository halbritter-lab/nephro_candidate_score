# Title: gnomad gene constraints

# load libraries
library(tidyverse)
library(R.utils)

# download and unzip gnomad constraint data
download_url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

download.file(url = download_url, 
              destfile = "gene_score/features/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")

gunzip(filename = "gene_score/features/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", 
       destname = "gene_score/features/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt")

# load and select gnomad gene constraints
gnomad_constraints <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene.txt") %>%
  dplyr::select(gene, 
                gene_id,
                obs_mis,
                exp_mis,
                oe_mis,
                mu_mis,
                possible_mis,
                obs_mis_pphen,
                exp_mis_pphen,
                oe_mis_pphen,
                possible_mis_pphen,
                obs_syn,
                exp_syn,
                oe_syn,
                mu_syn,
                possible_syn,
                obs_lof,
                mu_lof,
                possible_lof,
                exp_lof,
                pLI,
                pRec,
                pNull,
                oe_lof,
                syn_z, 
                mis_z, 
                lof_z,
                oe_lof_upper_rank,
                n_sites,
                classic_caf,
                max_af,
                p,
                exp_hom_lof,
                cds_length,
                num_coding_exons,
                gene_length
                )

# add prefix
names(gnomad_constraints)[3:length(names(gnomad_constraints))] <- paste0("gnomad_", names(gnomad_constraints)[3:length(names(gnomad_constraints))])

# write results
write.csv(gnomad_constraints, paste0("gene_score/features/results/gnomad_constraints_", creation_date, ".csv"), row.names = FALSE)
