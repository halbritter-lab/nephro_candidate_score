# Title: Average observed to expected CpG ratio of exons region of canonical transcripts

# load libraries
library(tidyverse)
library(biomaRt)
library(progress)
library(BSgenome.Hsapiens.UCSC.hg38) # approx 700MB
library(config)
library(R.utils)

# define relative script path
project_topic <- "nephrology"
project_name <- "nephro_candidate_score"
script_path <- "/gene_score/features/"

# read configs
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"),
                           config = project_topic)

# set working directory
setwd(paste0(config_vars$projectsdir, project_name, script_path))

# load canonical transcripts (from script "promoter_CpG_o2e_ratio.R")
canon_ts <- read_csv(paste0("results/ensembl_canonical_ts_", config_vars$creation_date, ".csv.gz"),
                             show_col_types = FALSE,
                             na = c("NA", "NaN", " ", ""))

# download Ensembl data 
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      version = config_vars$ensembl_biomart_version)

exons_coordinates <- getBM(attributes = c("ensembl_gene_id",
                                          "external_gene_name",
                                          "ensembl_transcript_id",
                                          "chromosome_name",
                                          "exon_chrom_start",
                                          "exon_chrom_end",
                                          "cds_start",
                                          "cds_end",
                                          "genomic_coding_start",
                                          "genomic_coding_end"),
                           filters = "ensembl_transcript_id",
                           values = list(canon_ts$ensembl_transcript_id),
                           mart = ensembl) %>% 
  filter(!is.na(cds_start)) %>% # filter out non-coding exons
  mutate(all_coding = case_when((cds_end - cds_start) == (exon_chrom_end - exon_chrom_start) ~ TRUE, TRUE ~ FALSE)) %>% 
  dplyr::rename(symbol = external_gene_name)


# create a GRanges object
canonical_exons_granges <- GRanges(
  seqnames = paste0("chr", exons_coordinates$chromosome_name),
  ranges = IRanges(start = exons_coordinates$genomic_coding_start, end = exons_coordinates$genomic_coding_end)
)

# retrieve sequences from the genome
exons_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, canonical_exons_granges) %>% as.character()

# join with former df
exons_seq <- cbind(exons_coordinates, data.frame(sequence = exons_seq))

# calculate CpG observed to expected ratio of each exon
calculate_CpG_obs_to_exp_ratio <- function(seq){
  seq_length <- str_length(seq)
  CG <- str_count(seq, "CG")
  C <- str_count(seq, "C")
  G <- str_count(seq, "G")
  CpG_obs_to_exp_ratio <- (CG / (C*G))*seq_length
  return(CpG_obs_to_exp_ratio)
}

exons_ratios <- exons_seq %>% 
  mutate(CpG_o2e_ratio = calculate_CpG_obs_to_exp_ratio(sequence))

# calculate average CpG-o2e-ratio per gene
exons_ratios <- exons_ratios %>% 
  group_by(ensembl_gene_id, symbol) %>% 
  summarise(avg_exon_CpG_o2e_ratio = mean(CpG_o2e_ratio, na.rm = TRUE))


# write results
write.csv(exons_ratios, 
          paste0("results/canonical_ts_exons_CpG_obs_to_exp_ratio_" , config_vars$creation_date, ".csv"), 
          row.names=FALSE)

gzip(paste0("results/canonical_ts_exons_CpG_obs_to_exp_ratio_" , config_vars$creation_date, ".csv"),
     overwrite = TRUE)
