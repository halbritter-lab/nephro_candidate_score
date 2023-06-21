# CAVE: unfinished => check TODO's !!

# Title: Exon and promoter conservation scores and coding sequence lengths

# load libraries
library(tidyverse)
library(biomaRt)
library(progress)
library(phastCons100way.UCSC.hg38)  # PhastCons score ranges from 0 to 1 and represents the probability that a given nucleotide is conserved

# load canonical transcripts (from script "promoter_CpG_o2e_ratio.R")
canon_ts <- read.csv(paste0("gene_score/features/results/ensembl_canonical_ts_", creation_date, ".csv"),
                       na.strings = c("NA", "NaN", " ", ""))

# download Ensembl data 
bm_version <- 109  # Ensembl Biomart Genes version
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      version = bm_version)

exons_coordinates <- getBM(attributes = c("ensembl_gene_id",
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
  filter(!is.na(cds_start)) %>% # filter non-coding exons
  mutate(all_coding = case_when((cds_end - cds_start) == (exon_chrom_end - exon_chrom_start) ~ TRUE, TRUE ~ FALSE)) 

## Exon conservation scores
# function to return average of phastCons scores of a given transcript
# loads phastCons scores per transcript and not per exon => faster
avg_phastCons_score_per_transcript <- function(transcript_id){
  trans_df <- exons_coordinates %>%
    filter(ensembl_transcript_id == transcript_id)
  scores <- gscores(phastCons100way.UCSC.hg38, 
                    GRanges(seqnames = paste0("chr", trans_df$chromosome_name), 
                             IRanges(start = trans_df$genomic_coding_start, end = trans_df$genomic_coding_end))) %>%
    data.frame() %>% 
    mutate(sum_phastCon_per_exon = default * width)  # calculate exon sum of phastCon scores back from exon average
  
  # calculate average phastCons score per transcript
  avg <- sum(scores$sum_phastCon_per_exon)/sum(scores$width)
             
  return(avg)
}

# get average phasCons scores for each canonical transcript (long calculation time)
pb <- progress_bar$new(total = length(unique(exons_coordinates$ensembl_gene_id)))

avg_phasCons_scores <- exons_coordinates %>%  
  dplyr::select(ensembl_gene_id, ensembl_transcript_id) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(avg_phasCons = {
    pb$tick()
    avg_phastCons_score_per_transcript(ensembl_transcript_id)
    })

# write results
write.csv(avg_phasCons_scores, 
          paste0("gene_score/features/results/avg_phasCons_scores_per_transcript_" , creation_date, ".csv"), 
          row.names=FALSE)



## Promoter conservation scores
# function to calculate avearge phastCons score of promoter
avg_phastCons_score_promoter <- function(chrom, start, end){
   scores <-  gscores(phastCons100way.UCSC.hg38, GRanges(seqnames=paste0("chr", chrom), IRanges(start=start:end, width=1))) %>%
    data.frame() %>% 
     .$default %>% 
     mean(na.rm=TRUE)
  return(scores)
}

# get average phasCons scores for each promoter (long calculation time)
pb <- progress_bar$new(total = length(unique(canon_trans$ensembl_gene_id)))
avg_phasCons_prom <- canon_trans %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(avg_phasCons_promoter = {
    pb$tick()
    avg_phastCons_score_promoter(chrom = chrom, start = start, end = end)
  })

# write results
write.csv(avg_phasCons_prom, 
          paste0("gene_score/features/results/avg_phasCons_promoter_" , creation_date, ".csv"), 
          row.names=FALSE)


