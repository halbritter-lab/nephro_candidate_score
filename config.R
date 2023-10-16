# configuration file

PROJECT_DIR <- "/Users/nrank/Desktop/BioInf/halbritter/nephro_candidate_score/" # set your project directory here
creation_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

kidney_genetics_version <- "2023-10-04"  # date of kidney-genetics version, check on https://github.com/halbritter-lab/kidney-genetics/tree/main/analyses/merged for latest version  

hgnc_gt_version <- "2023-06-21" # date of most recent HGNC annotated gene table version, https://github.com/halbritter-lab/kidney-genetics/tree/main/analyses/A_AnnotationHGNC/results 


omim_download_url <- "https://data.omim.org/downloads/1R83aw9QTZ66GmeYgYPVAQ/genemap2.txt" # set your omim download link for genemap


# setting up directories etc.
setwd(PROJECT_DIR)
