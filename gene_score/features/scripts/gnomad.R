# Title: gnomad gene constraints

# load libraries
library(tidyverse)
library(utils)
library(R.utils)
library(config)


# define relative script path
project_topic <- "nephrology"
project_name <- "nephro_candidate_score"
script_path <- "/gene_score/features/"

# read configs
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"),
                           config = project_topic)

# set working directory
setwd(paste0(config_vars$projectsdir, project_name, script_path))

# download and unzip gnomad constraint data
download_url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

download.file(url = download_url, 
              destfile = "raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")

gunzip(filename = "raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", 
       destname = "raw/gnomad.v2.1.1.lof_metrics.by_gene.txt")

# load and select gnomad gene constraints
gnomad_constraints <- read.delim("raw/gnomad.v2.1.1.lof_metrics.by_gene.txt") %>%
  dplyr::select(ensembl_gene_id = gene_id, 
                symbol = gene,
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
write.csv(gnomad_constraints, 
          paste0("results/gnomad_constraints_", config_vars$creation_date_gs, ".csv"), 
          row.names = FALSE)

gzip(paste0("results/gnomad_constraints_", config_vars$creation_date_gs, ".csv"),
     overwrite = TRUE)
