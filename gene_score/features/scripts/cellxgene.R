# Title: ssRNA data from cellxgene

# load libraries
library(tidyverse)
library(httr)
library(jsonlite)

# set the URL
query_url <- "https://api.cellxgene.cziscience.com/wmg/v1/query"

# set the header
headers <- c("Content-Type" = "application/json")

# set the request body - NOTE: following code should be adapted if multiple organisms are queried at once
body <- list(
  filter = list(
    dataset_ids = c("0b4a15a7-4e9e-4555-9733-2423e5c66469"),
    development_stage_ontology_term_ids = list(),
    disease_ontology_term_ids = list(),
    gene_ontology_term_ids = c("ENSG00000089597", "ENSG00000000005"),
    organism_ontology_term_id = unbox("NCBITaxon:9606"),
    self_reported_ethnicity_ontology_term_ids = list(),
    sex_ontology_term_ids = list(),
    tissue_ontology_term_ids = c("UBERON:0002113")
  ),
  is_rollup = unbox(TRUE)
)

# Send the POST request
response <- POST(query_url, add_headers(headers), body = toJSON(body))

# Get the response content
content <- content(response, as="text")

# transform json content into a data frame with genes as rows and expression values as columns
expr_vals <- fromJSON(content)

tissue <- "UBERON:0002113" # TODO change

expr_lists <- sapply(expr_vals[["expression_summary"]], "[[", tissue) # one list with lists of expression values, elements correspond to genes

# combine lists of expression values into a data frame
expr_df <- dplyr::bind_rows(lapply(expr_lists, as.data.frame))
gene_id <- names(expr_vals[["expression_summary"]]) # list of genes in the json file
expr_df <- cbind(gene_id, expr_df)

# TODO: maybe remove "aggregated" from column names

View(expr_df)

# CONTINUE HERE: extract all genes from primary filter

#############################################
# Extract primary filter dimensions
# set the URL
prim_filt_url <- "https://api.cellxgene.cziscience.com/wmg/v1/primary_filter_dimensions"

# set the header
headers <- c("Content-Type" = "application/json")

# send the GET request
response <- GET(prim_filt_url, add_headers(headers))

# get the response content
content <- content(response, as="text")

# x <- fromJSON(content)
x <- read_json(content)

# https://www.geeksforgeeks.org/how-to-read-large-json-file-in-r/

setwd("/Users/nrank/Desktop/BioInf/halbritter/nephro_gene_score/")

writeLines(content[[1]], "junk_content.txt")                 # Apply writeLines function

data <- jsonlite::read_json("junk_content.txt")

head(data, 3)



