# Title:  Observed to expected CpG ratio of promoter region of canonical transcripts

# load libraries
library(biomaRt) # http://www.ensembl.org/info/data/biomart/biomart_r_package.html#biomartexamples
library(tidyverse)
library(jsonlite)
library(BSgenome.Hsapiens.UCSC.hg38) # approx 700MB

# download Ensembl data for all protein coding transcripts in Homo sapiens, Human genes (GRCh38.p13), GRCh38.p13

# download Ensembl data 
bm_version <- 109  # Ensembl Biomart Genes version
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      version = bm_version)

mart_exp <- getBM(attributes = c('ensembl_gene_id',
                                 'ensembl_transcript_id',
                                 'ensembl_transcript_id_version',
                                 'hgnc_id',
                                 'chromosome_name',
                                 'transcript_start',
                                 'transcript_end',
                                 'transcription_start_site',
                                 'transcript_is_canonical',
                                 'transcript_biotype'),
                  filter = 'transcript_biotype', 
                  values = 'protein_coding',
                  mart = ensembl) %>% 
  filter(chromosome_name %in% c(as.character(seq(1:22)), "X", "Y"))

# Coordinates promoters canonical transcripts: 4kb regions surrounding the TSS (+-2kb)
canonical_prom <- mart_exp %>% 
  filter(transcript_is_canonical == 1) %>% 
  mutate(start = transcription_start_site - 2000,
         end = transcription_start_site + 2000) %>% 
  dplyr::select(ensembl_gene_id, ensembl_transcript_id, chrom=chromosome_name, start, end) %>% 
  distinct()

# write canonical transcripts
write_csv(canonical_prom, 
          paste0("gene_score/features/results/ensembl_canonical_ts_" , creation_date, ".csv"))

# create a GRanges object
canonical_prom_granges <- GRanges(
  seqnames = paste0("chr", canonical_prom$chrom),
  ranges = IRanges(start = canonical_prom$start, end = canonical_prom$end)
)
# retrieve sequences from the genome
prom_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, canonical_prom_granges) %>% 
  as.character()

# join with former df
canonical_prom <- cbind(canonical_prom, data.frame(sequence = prom_seq))

# calculate CpG observed to expected ratio 
calculate_CpG_obs_to_exp_ratio <- function(seq){
  seq_length <- str_length(seq)
  CG <- str_count(seq, "CG")
  C <- str_count(seq, "C")
  G <- str_count(seq, "G")
  CpG_obs_to_exp_ratio <- (CG / (C*G))*seq_length
  return(CpG_obs_to_exp_ratio)
}

canonical_prom <- canonical_prom %>% 
  mutate(CpG_o2e_ratio = calculate_CpG_obs_to_exp_ratio(sequence))

# write results
write.csv(canonical_prom[, c("ensembl_gene_id", "CpG_o2e_ratio")], 
          paste0("gene_score/features/results/canonical_promoter_CpG_obs_to_exp_ratio_" , creation_date, ".csv"), 
          row.names=FALSE)

