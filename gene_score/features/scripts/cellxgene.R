# Title: ssRNA data from cellxgene

# load libraries
library(tidyverse)
library(httr)
library(jsonlite)
library(progress)

########## DATASETS ##################
# get all datasets from cellxgene
datasets_url <- "https://api.cellxgene.cziscience.com/dp/v1/datasets/index"

# set the header
headers <- c("Content-Type" = "application/json")

# get the response
response <- GET(datasets_url, add_headers(headers))

# get the response content
content <- content(response, as = "text", encoding = "UTF-8")

# write json file
writeLines(content[[1]], paste0("gene_score/features/raw/cellxgene_datasets_", creation_date, ".json"))

# read json file
datasets_list <- jsonlite::read_json(paste0("gene_score/features/raw/cellxgene_datasets_", creation_date, ".json"))

# wrap list in a tibble
datasets <- tibble(ds = datasets_list)

# create a summary with most important info of each dataset
# TODO: maybe collapse organisms and diseases?
datasets_df <- datasets %>%
  unnest_wider(ds) %>%
  dplyr::select(id, collection_id, name, organism, cell_count, disease, published, explorer_url) %>% 
  unnest_longer(organism) %>%  
  hoist(organism, 
        organism_label = "label", 
        organism_ontology_term_id = "ontology_term_id"
        ) %>% 
  unnest_longer(disease) %>% 
  hoist(disease,
        disease_label = "label",
        disease_ontology_term_id = "ontology_term_id"
        ) %>% 
  mutate(query_id = regmatches(explorer_url, regexpr("[[:alnum:]]{8}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}", explorer_url))
)

# write results
write.csv(datasets_df, paste0("gene_score/features/results/cellxgene_datasets_summary", creation_date, ".csv"), row.names = FALSE)


########## PRIMARY FILTER DIMENSIONS ##################
# set the URL
prim_filt_url <- "https://api.cellxgene.cziscience.com/wmg/v1/primary_filter_dimensions"

# set the header
headers <- c("Content-Type" = "application/json")

# send the GET request
response <- GET(prim_filt_url, add_headers(headers))

# get the response content
content <- content(response, as = "text", encoding = "UTF-8")

# save content (workaround as fromJSON() is too slow, read_json() reads the file a memory-efficient way)
writeLines(content[[1]], paste0("gene_score/features/raw/primary_filter_dim_content_", creation_date, ".json"))            

# read in json file line by line
prim_filt_dim_list <- jsonlite::read_json(paste0("gene_score/features/raw/primary_filter_dim_content_", creation_date, ".json"))

# wrap list in a tibble
prim_filt_dim <- tibble(pfd = prim_filt_dim_list)

# primary filter dimensions: organism_terms
pfd_ot <- prim_filt_dim_list$organism_terms  

pfd_ot_df <- imap_dfr(pfd_ot, ~ tibble(taxon_id = names(.x), organism = unlist(.x)))

# primary filter dimensions: gene_terms
pfd_gt <- prim_filt_dim_list$gene_terms 

pfd_gt_df <- pfd_gt %>%
  imap_dfr(~ map_df(.x, ~ tibble(ensembl_id = names(.x), symbol = as.character(.x)))) %>%
  mutate(taxon_id = rep(names(pfd_gt), lengths(pfd_gt)))

# primary filter dimensions: tissue_terms
pfd_tt <- prim_filt_dim_list$tissue_terms 

pfd_tt_df <- pfd_tt %>%
  imap_dfr(~ map_df(.x, ~ tibble(uberon_id = names(.x), tissue = as.character(.x)))) %>%
  mutate(taxon_id = rep(names(pfd_tt), lengths(pfd_tt)))



########## QUERY EXPRESSION VALUES ##########
# TODO: customize, rewrite etc. 

# set the URL
query_url <- "https://api.cellxgene.cziscience.com/wmg/v1/query"

# set the header
headers <- c("Content-Type" = "application/json")

# functions for querying expression values from cellxgene
# TODO: if needed, customize for more than one tissue etc.
create_query_body <- function(datasets_vec,
                              query_taxon,
                              query_tissue,
                              development_stage_list = list(),
                              disease_list = list(),
                              gene_list = list(),
                              self_reported_ethnicity_list = list(),
                              sex_list = list()
) {
  body <- list(
    filter = list(
      dataset_ids = datasets_vec,
      development_stage_ontology_term_ids = development_stage_list,
      disease_ontology_term_ids = disease_list,
      gene_ontology_term_ids = gene_list,
      organism_ontology_term_id = unbox(query_taxon),
      self_reported_ethnicity_ontology_term_ids = self_reported_ethnicity_list,
      sex_ontology_term_ids = sex_list,
      tissue_ontology_term_ids = c(query_tissue)
    ),
    is_rollup = unbox(TRUE) #TODO: change?
  )
  return(body)
}

query_cellxgene_POST <- function(url, body) {
  # send the POST request
  response <- POST(url, add_headers(headers), body = toJSON(body))
  
  # set the response content
  content <- content(response, as = "text", encoding = "UTF-8")

  # parse the content
  content_list <- fromJSON(content)
  
  return(content_list)
}


list_to_expr_df <- function(content_list) {
  expr_vals_df <- imap_dfr(content_list$expression_summary, ~ tibble(
    ensembl_id = .y,
    uberon_id = names(.x),
    cell_id = names(.x[[1]]),
    me = unname(unlist(lapply(lapply(.x[[1]], "[[", 1), "[[", "me"))),
    n = unname(unlist(lapply(lapply(.x[[1]], "[[", 1), "[[", "n"))),
    pc = unname(unlist(lapply(lapply(.x[[1]], "[[", 1), "[[", "pc"))),
    tpc = unname(unlist(lapply(lapply(.x[[1]], "[[", 1), "[[", "tpc")))
  )) 
  
  # pivot df
  expr_vals_piv <- expr_vals_df %>%
    pivot_wider(
      names_from = cell_id,
      values_from = c("me", "n", "pc", "tpc"),
      names_glue = "{cell_id}_{.value}"
    )
  return(expr_vals_piv)
}


get_cellxgene_expr_val <- function(datasets_vec,
                                   query_taxon,
                                   query_tissue,
                                   development_stage_list = list(),
                                   disease_list = list(),
                                   gene_ontology_term_ids = list(),
                                   self_reported_ethnicity_ontology_term_ids = list(),
                                   sex_ontology_term_ids = list(),
                                   url
                                   ) {
  body <- create_query_body(datasets_vec,
                            query_taxon,
                            query_tissue,
                            development_stage_list,
                            disease_list,
                            gene_ontology_term_ids,
                            self_reported_ethnicity_ontology_term_ids,
                            sex_ontology_term_ids)
  
  content_list <- query_cellxgene_POST(url, body)
  
  if (length(content_list$expression_summary) >= 1) {
    res_df <- list_to_expr_df(content_list)
    return(res_df)
  }else {
    return(data.frame())
    }
}
                                   


###### GET EXPRESSION VALUES FOR "0b4a15a7-4e9e-4555-9733-2423e5c66469" ######
ds_id <- c("0b4a15a7-4e9e-4555-9733-2423e5c66469")
query_taxon <- "NCBITaxon:9606"
query_tissue <- "UBERON:0002113"
gene_vec <- pfd_gt_df %>%
  filter(taxon_id %in% query_taxon, ensembl_id %in% HGNC_table$ensembl_gene_id) %>%
  .$ensembl_id %>%
  as.vector()

# NOTE: querying all genes at once is not possible => split up the gene vector into junks of approximately length 1000
gene_split_list <- split(gene_vec, ceiling(seq_along(gene_vec) / 1000))


# get expression values chunk wise
expr_val_combined <- data.frame()
pb <- progress_bar$new(total = length(gene_split_list))


for (i in gene_split_list){
  pb$tick()  # Increment progress bar
  cat("\n")
  sub_df <- get_cellxgene_expr_val(datasets_vec = ds_id,
                                 query_taxon = query_taxon,
                                 query_tissue = query_tissue,
                                 development_stage_list = list(),
                                 disease_list = list(),
                                 gene_ontology_term_ids = i,
                                 self_reported_ethnicity_ontology_term_ids = list(),
                                 sex_ontology_term_ids = list(),
                                 url = "https://api.cellxgene.cziscience.com/wmg/v1/query")
  
  expr_val_combined <- dplyr::bind_rows(expr_val_combined, sub_df)
}

expr_val_combined <- expr_val_combined %>%
  dplyr::select(ensembl_id, ends_with("_me"), ends_with("_pc"))

# write results
write.csv(expr_val_combined, paste0("gene_score/features/results/cellxgene_expr_0b4a15a7-4e9e-4555-9733-2423e5c66469_", creation_date, ".csv"), row.names = FALSE)

###### GET EXPRESSION VALUES FOR "d7dcfd8f-2ee7-4385-b9ac-e074c23ed190" (FETAL KIDNEY) ######
ds_id <- c("d7dcfd8f-2ee7-4385-b9ac-e074c23ed190")
query_taxon <- "NCBITaxon:9606"
query_tissue <- "UBERON:0002113"
gene_vec <- pfd_gt_df %>%
  filter(taxon_id %in% query_taxon, ensembl_id %in% HGNC_table$ensembl_gene_id) %>%
  .$ensembl_id %>%
  as.vector()

# NOTE: querying all genes at once is not possible => split up the gene vector into junks of approximately length 1000
gene_split_list <- split(gene_vec, ceiling(seq_along(gene_vec) / 1000))

# get expression values chunk wise
expr_val_combined <- data.frame()
pb <- progress_bar$new(total = length(gene_split_list))

for (i in gene_split_list){
  pb$tick()  # Increment progress bar
  cat("\n")
  sub_df <- get_cellxgene_expr_val(datasets_vec = ds_id,
                                   query_taxon = query_taxon,
                                   query_tissue = query_tissue,
                                   development_stage_list = list(),
                                   disease_list = list(),
                                   gene_ontology_term_ids = i,
                                   self_reported_ethnicity_ontology_term_ids = list(),
                                   sex_ontology_term_ids = list(),
                                   url = "https://api.cellxgene.cziscience.com/wmg/v1/query")
  
  expr_val_combined <- dplyr::bind_rows(expr_val_combined, sub_df)
}

expr_val_combined <- expr_val_combined %>%
  dplyr::select(ensembl_id, ends_with("_me"), ends_with("_pc"))

# write results
write.csv(expr_val_combined, paste0("gene_score/features/results/cellxgene_expr_d7dcfd8f-2ee7-4385-b9ac-e074c23ed190_", creation_date, ".csv"), row.names = FALSE)





