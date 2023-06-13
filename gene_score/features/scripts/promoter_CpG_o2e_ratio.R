# Title: Expected to observed CpG ratio of promoter region of canonical transcripts

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

ensembl_transcripts <- getBM(attributes=c('ensembl_gene_id',
                                     'ensembl_transcript_id',
                                     'ensembl_transcript_id_version',
                                     'hgnc_id',
                                     'hgnc_symbol',
                                     'chromosome_name',
                                     'start_position',
                                     'end_position',
                                     'transcript_start',
                                     'transcript_end',
                                     'transcription_start_site',
                                     'transcript_gencode_basic',
                                     'transcript_appris',
                                     'transcript_mane_select',
                                     'transcript_biotype'),
                        filter='transcript_biotype', 
                        values = 'protein_coding',
                        mart=ensembl) %>% 
  filter(chromosome_name %in% c(as.character(seq(1:22)), "X", "Y")) %>% 
  mutate(length = transcript_end - transcript_start)

# TODO: replace HGNC_OMIM from Leitao with our table, also further down!!!
HGNC_OMIM <- read_csv("~/Desktop/BioInf/halbritter/nephro_gene_score/Leitao/christopher-schroeder-chrX_gene_predictions-1ebb650/additional/gene_annotations/results/HGNC_OMIM_all_chr/HGNCgenes_OMIMannotated.csv") %>% 
  rename(ensemblID = `Ensembl gene ID`)

## MANE and APPRIS annotations,  #TODO: add citation Leitao
# - MANE select transcripts (Matched Annotation between NCBI and EBI) were independently identified by both Ensembl and NCBI as the most biologically relevant.   
# - APPRIS is a system to annotate alternatively spliced transcripts based on a range of computational methods.

MANE_transcripts <- ensembl_transcripts %>% 
  filter(!is.na(transcript_mane_select) & ensembl_gene_id %in% HGNC_OMIM$ensemblID)

# TODO: from APPROS use principal 1, principal 2 alternative etc? or only lenghts? (Here only length)
APPRIS_transcripts <-  ensembl_transcripts %>% 
  filter(!(ensembl_gene_id %in% MANE_transcripts$ensembl_gene_id) &
           !is.na(transcript_appris) &
           ensembl_gene_id %in% HGNC_OMIM$ensemblID) %>%
  group_by(ensembl_gene_id) %>%
  arrange(desc(length)) %>%
  dplyr::slice(1)

other_transcripts <-  ensembl_transcripts %>%
  filter(ensembl_gene_id %in% HGNC_OMIM$ensemblID &
           !(ensembl_gene_id %in% MANE_transcripts$ensembl_gene_id) &
           !(ensembl_gene_id) %in% APPRIS_transcripts$ensembl_gene_id) %>%
  group_by(ensembl_gene_id) %>%
  arrange(desc(length)) %>%
  dplyr::slice(1)


# get canonical transcripts (Leitao et al.) according to the following order:
# 1 - MANE transcript if it exists (vast majority. approx. 80% genes)  
# 2 - longest transcript with APPRIS annotation (approx. 20% genes)  
# 3 - existing transcript (here none)  

canonical_transcripts <-  MANE_transcripts %>%
  full_join(APPRIS_transcripts) %>%
  full_join(other_transcripts) 

# get promoter coordinates of canonical transcripts: 4kb regions surrounding the TSS (+-2kb)
canonical_TSS <-  canonical_transcripts %>%
  mutate(start = transcription_start_site - 2000,
         end = transcription_start_site + 2000) %>%
  dplyr::select(ensembl_gene_id, canon_transcript_id = ensembl_transcript_id, chrom=chromosome_name, start, end ) 

# create a GRanges object
canonical_TSS_granges <- GRanges(
  seqnames = paste0("chr", canonical_TSS$chrom),
  ranges = IRanges(start = canonical_TSS$start, end = canonical_TSS$end)
)

# retrieve sequences from BSgenome.Hsapiens.UCSC.hg38
prom_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, canonical_TSS_granges) %>% 
  as.character()

# join with former df
canonical_TSS <- cbind(canonical_TSS, data.frame(sequence = prom_seq))

# calculate CpG observed to expected ratio 
calculate_CpG_obs_to_exp_ratio <- function(seq){
  seq_length <- str_length(seq)
  CG <- str_count(seq, "CG")
  C <- str_count(seq, "C")
  G <- str_count(seq, "G")
  CpG_obs_to_exp_ratio <- (CG / (C*G))*seq_length
  return(CpG_obs_to_exp_ratio)
}

canonical_TSS <- canonical_TSS %>% 
  mutate(CpG_o2e_ratio = calculate_CpG_obs_to_exp_ratio(sequence))

# write results
write.csv(canonical_TSS[, c("ensembl_gene_id", "CpG_o2e_ratio")], 
          paste0("gene_score/features/results/canonical_promoter_CpG_obs_to_exp_ratio_" , creation_date, ".csv"), 
          row.names=FALSE)

