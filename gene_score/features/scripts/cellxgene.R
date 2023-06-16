# Title: ssRNA data from cellxgene

# Title: ssRNA data from cellxgene

# load libraries
library(tidyverse)
library(httr)
library(jsonlite)

setwd("/Users/nrank/Desktop/BioInf/halbritter/nephro_gene_score/cellxgene")
creation_date <- "2023-04-21"


########## DATASETS ##################
# get all datasets from cellxgene
datasets_url <- "https://api.cellxgene.cziscience.com/dp/v1/datasets/index"

# set the header
headers <- c("Content-Type" = "application/json")

# get the response
response <- GET(datasets_url, add_headers(headers))

# get the response content
content <- content(response, as="text")

# content

# write json file
writeLines(content[[1]], paste0("cellxgene_datasets_", creation_date, ".json"))   

# read json file
datasets <- jsonlite::read_json(paste0("cellxgene_datasets_", creation_date, ".json"), flatten=T) #try this
datasets <- jsonlite::read_json(paste0("cellxgene_datasets_", creation_date, ".json")) 

# https://medium.com/@Periscopic/cozy-collecting-part-2-5e717588e37b # try this!!!

# check how many organisms a dataset has
sapply(sapply(datasets, "[[", "organism"), length)

datasets[[583]][['organism']] # has three organism, "367d95c0-0eb0-4dae-8276-9407239421ee"
sapply(datasets[[583]][["organism"]], "[[", "label")
d583 <- datasets[[583]]

d583$cell_count
d583$assay

# check how many cellcounts a dataset has
sapply(sapply(datasets, "[[", "cell_count"), length) %>% unique() # only 1

# check how many collection_id a dataset has
sapply(sapply(datasets, "[[", "collection_id"), length) %>% unique() # only 1

# check multivalue fields in datasets
for (i in names(d583)){
  print(i)
  a <- sapply(sapply(datasets, "[[", i), length) %>% unique() 
  print(a)
}
# multivalue fields are: "assay", "cell_type", "cell_type_ancestors", "dataset_assets", 
# "development_stage", "development_stage_ancestors", "disease", "donor_id", "organism", "self_reported_ethnicity", "sex",
# "suspension_type", "tissue", "tissue_ancestors"

id <- sapply(datasets, "[[", "id")
cell_count <- sapply(datasets, "[[", "cell_count")
collection_id <- sapply(datasets, "[[", "collection_id")

# CONTINUE HERE 
organisms <- sapply(sapply(datasets, "[[", "organism"), "[[", "label")

set3 <- organisms[c(1,2,583)]


sapply(d583$organism, "[[", "ontology_term_id")

index_function <- function(x){
  
  
})



org <- c()
org_length <- c()
for (i in seq(length(datasets))){
  print(i)
  ol <- length(datasets[[i]]["organsim"])
  orga <- datasets[[i]][['organism']][[1]][['label']]
  print(orga)
  org <- c(org, orga)
  org_length <- c(org_length, ol)
}

org %>% as.data.frame() %>% rename(organism = ".") %>% View
org %>% as.data.frame() %>% rename(organism = ".") %>% group_by(organism) %>%  summarize(count=n()) %>% View




datasets$filter_dims$datasets %>% as.data.frame() %>% View

datasets$filter_dims %>% View



names(datasets$filter_dims)
# [1] "cell_type_terms"               "datasets"                      "development_stage_terms"       "disease_terms"                
# [5] "self_reported_ethnicity_terms" "sex_terms"                     "tissue_terms"      

# names(datasets$filter_dims$datasets)
# [1] "collection_id"    "collection_label" "id"               "label" 







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
# Extract filter dimensions
# set the URL
filt_dim_url <- "https://api.cellxgene.cziscience.com/wmg/v1/filters"

# set the header
headers <- c("Content-Type" = "application/json")


# set the request body - NOTE: following code should be adapted if multiple organisms are queried at once
body <- list(
  filter = list(
    dataset_ids = list(),
    development_stage_ontology_term_ids = list(),
    disease_ontology_term_ids = list(),
    gene_ontology_term_ids = list(),
    # organism_ontology_term_id = unbox("NCBITaxon:9606&NCBITaxon.10090"),
    
    organism_ontology_term_id = unbox("NCBITaxon.10090"),
    self_reported_ethnicity_ontology_term_ids = list(),
    sex_ontology_term_ids = list(),
    tissue_ontology_term_ids = list()
  ),
  is_rollup = unbox(TRUE)
)

# Send the POST request
response <- POST(filt_dim_url, add_headers(headers), body = toJSON(body))

# get the response content
content <- content(response, as="text")

content
writeLines(content[[1]], paste0("junk_content3.json"))            

data <- fromJSON(content)
data$filter_dims$datasets$collection_id %>% length





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

# https://www.geeksforgeeks.org/how-to-read-large-json-file-in-r/

setwd("/Users/nrank/Desktop/BioInf/halbritter/nephro_gene_score/")
creation_date <- "2023-04-21"

# save content (workaround as fromJSON() is too slow, read_json() reads the file a memory-efficient way)
writeLines(content[[1]], paste0("primary_filter_dim_content_", creation_date, ".json"))            

# read in json file line by line
prim_filt_dim <- jsonlite::read_json(paste0("primary_filter_dim_content_", creation_date, ".json"))

# function to return organism dimensions
get_organism_dim <- function(prim_filt_dim_df){
  organism_dim <- prim_filt_dim$organism_terms %>% 
    data.frame() %>% 
    t() %>%  # Value = matrix
    as.data.frame() %>% 
    rownames_to_column(var = "taxon_id") %>% 
    rename(organsim = V1)
  return(organism_dim)
}

# get organism dimensions
organism_dim <- get_organism_dim(prim_filt_dim)

# function to return gene dimensions for a given taxon
get_gene_dim <- function(prim_filt_dim_df, taxon_id){
  gene_dim <- prim_filt_dim_df$gene_terms[[taxon_id]] %>% 
    data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "ensembl_id") %>% 
    rename(symbol = V1)
  return(gene_dim)
}

# get gene dimensions for Homo sapiens (NCBITaxon:9606) 
NCBITaxon_9606_gene_dim <- get_gene_dim(prim_filt_dim, "NCBITaxon:9606")

#


gene_dim <- prim_filt_dim$gene_terms[[taxon]] %>% 
  data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ensembl_id") %>% 
  rename(symbol = V1)




# names_df <- names(junk$expression_summary)
# 
# df1 <- data.frame(a = c(1:5), b = c(6:10))
# df2 <- data.frame(a = c(11:15), b = c(16:20), c = LETTERS[1:5])
# df3 <- data.frame(a = c(10:12), b = c(6:8))
# 
# dplyr::bind_rows(df1, df2, df3)
# 
# 
# 
# sapply(oneg[[tissue]], "[[", "aggregated")
# 
# 
# lapply(oneg, .[tissue])
#   
#   
# junk$expression_summary$ENSG00000089597$`UBERON:0002113`[c("CL:0000066", "CL:0002306")]
# 
# b <- junk$expression_summary$ENSG00000089597$`UBERON:0002113`
# c <- b[c("CL:0000066", "CL:0000084", "CL:0000115")]
# 
# 
# junk$expression_summary$ENSG00000089597[tissue]
#   

# # works!
# body2 <- '{"filter":{"dataset_ids":["0b4a15a7-4e9e-4555-9733-2423e5c66469"],"development_stage_ontology_term_ids":[],"disease_ontology_term_ids":[],"gene_ontology_term_ids":["ENSG00000089597"],"organism_ontology_term_id":"NCBITaxon:9606","self_reported_ethnicity_ontology_term_ids":[],"sex_ontology_term_ids":[],"tissue_ontology_term_ids":["UBERON:0002113"]},"is_rollup":true}'
# response <- POST(url, add_headers(headers), body = body2)
# 
# # Get the response content
# content <- content(response, as="text")
# 
# # Print the response
# print(content)
# 
# 
# # body3
# body3 <- list(
#   filter = list(
#     dataset_ids = c("0b4a15a7-4e9e-4555-9733-2423e5c66469"),
#     development_stage_ontology_term_ids = list(),
#     organism_ontology_term_id = unbox("NCBITaxon:9606"),
#     tissue_ontology_term_ids = c("UBERON:0002113")
#   ),
#   is_rollup = TRUE
# )
# 
# toJSON(body3, auto_unbox=F)
# 
# 
# library(httr)
# library(jsonlite)
# 
# # Set the URL
# url <- "https://api.cellxgene.cziscience.com/wmg/v1/query"
# 
# # Set the headers
# headers <- c("Content-Type" = "application/json")
# 
# # Set the request body
# body <- list(
#   filter = list(
#     dataset_ids = list("0b4a15a7-4e9e-4555-9733-2423e5c66469"),
#     development_stage_ontology_term_ids = list(),
#     disease_ontology_term_ids = list(),
#     gene_ontology_term_ids = c("ENSG00000089597"),
#     organism_ontology_term_id = 'NCBITaxon:9606',
#     self_reported_ethnicity_ontology_term_ids = list(),
#     sex_ontology_term_ids = list(),
#     tissue_ontology_term_ids = list("UBERON:0002113")
#   ),
#   is_rollup = TRUE
# )
# 
# # Convert the body to JSON
# body_json <- toJSON(body)
# 
# # Send the POST request
# response <- POST(url, add_headers(headers), body = body_json)
# 
# # Get the response content
# content <- content(response, "text")
# 
# # Print the response
# print(content)
# 
# 
# 
# 
# 
# 
# 
# 
# # curl -v -H "Content-Type: application/json" -X POST -d  
# # '{"filter":{"dataset_ids":["0b4a15a7-4e9e-4555-9733-2423e5c66469"],"development_stage_ontology_term_ids":[],"disease_ontology_term_ids":[],"gene_ontology_term_ids":["ENSG00000089597"],
# # "organism_ontology_term_id":"NCBITaxon:9606","self_reported_ethnicity_ontology_term_ids":[],"sex_ontology_term_ids":[],"tissue_ontology_term_ids":["UBERON:0002113"]},"is_rollup":true}' 
# # https://api.cellxgene.cziscience.com/wmg/v1/query
# # 
# # 
# 
# 
# library(httr)
# 
# url <- "https://api.cellxgene.cziscience.com/wmg/v1/query"
# library(httr)
# 
# url <- "https://api.cellxgene.cziscience.com/wmg/v1/query"
# 
# headers <- c(
#   "Content-Type" = "application/json"
# )
# 
# body <- list(
#   filter = list(
#     dataset_ids = "0b4a15a7-4e9e-4555-9733-2423e5c66469",
#     development_stage_ontology_term_ids = list(),
#     disease_ontology_term_ids = list(),
#     gene_ontology_term_ids = "ENSG00000089597",
#     organism_ontology_term_id = "NCBITaxon:9606",
#     self_reported_ethnicity_ontology_term_ids = list(),
#     sex_ontology_term_ids = list(),
#     tissue_ontology_term_ids = "UBERON:0002113"
#   ),
#   is_rollup = TRUE
# )
# 
# response <- POST(url, body = toJSON(body), encode = "json", verbose(), add_headers(.headers = headers))
# 
# content(response)
# 
# 
# ###
# page_requested <- as.character(page_requested)
# 
# # get nonce
# natera_nonce <- natera_renasight_get_nonce()
# 
# library(curl)
# input_url <-  "https://api.cellxgene.cziscience.com/wmg/v1/query"
# 
# # open curl handle
# h <- new_handle()
# handle_setopt(h, customrequest = "POST")
# 
# # set curl form
# 
# body <- '{"filter":{"dataset_ids":["0b4a15a7-4e9e-4555-9733-2423e5c66469"],"development_stage_ontology_term_ids":[],"disease_ontology_term_ids":[],"gene_ontology_term_ids":["ENSG00000000005"],"organism_ontology_term_id":"NCBITaxon:9606","self_reported_ethnicity_ontology_term_ids":[],"sex_ontology_term_ids":[],"tissue_ontology_term_ids":["UBERON:0002113"]},"is_rollup":true}'
# 
# Test<-POST("https://api.cellxgene.cziscience.com/wmg/v1/query", body = toJSON(body))
# 
# 
# handle_setform(h, .list = list(
#   dataset_ids = "0b4a15a7-4e9e-4555-9733-2423e5c66469",
#   development_stage_ontology_term_ids = list(),
#   disease_ontology_term_ids = list(),
#   gene_ontology_term_ids = "ENSG00000089597",
#   organism_ontology_term_id = "NCBITaxon:9606",
#   self_reported_ethnicity_ontology_term_ids = list(),
#   sex_ontology_term_ids = list(),
#   tissue_ontology_term_ids = "UBERON:0002113"
# ))
# 
# r <- curl_fetch_memory(input_url, h)
# 
# 
# 
# handle_setform(h, .list = list(
#   dataset_ids = "0b4a15a7-4e9e-4555-9733-2423e5c66469",
#   development_stage_ontology_term_ids = "",
#   disease_ontology_term_ids = "",
#   gene_ontology_term_ids = "ENSG00000089597",
#   organism_ontology_term_id = "NCBITaxon:9606",
#   self_reported_ethnicity_ontology_term_ids = "",
#   sex_ontology_term_ids = "",
#   tissue_ontology_term_ids = "UBERON:0002113"
#   ))
# 
# 
# handle_setform(h, .list = list(filter = c(
#   dataset_ids = "0b4a15a7-4e9e-4555-9733-2423e5c66469",
#   development_stage_ontology_term_ids = list(),
#   disease_ontology_term_ids = list(),
#   gene_ontology_term_ids = "ENSG00000089597",
#   organism_ontology_term_id = "NCBITaxon:9606",
#   self_reported_ethnicity_ontology_term_ids = list(),
#   sex_ontology_term_ids = list(),
#   tissue_ontology_term_ids = "UBERON:0002113"
# )))
# 
# 
# library(curl)
# 
# library(curl)
# 
# url <- "https://api.cellxgene.cziscience.com/wmg/v1/query"
# headers <- c(
#   "Content-Type: application/json"
# )
# body <- '{"filter":{"dataset_ids":["0b4a15a7-4e9e-4555-9733-2423e5c66469"],"development_stage_ontology_term_ids":[],"disease_ontology_term_ids":[],"gene_ontology_term_ids":["ENSG00000089597"],"organism_ontology_term_id":"NCBITaxon:9606","self_reported_ethnicity_ontology_term_ids":[],"sex_ontology_term_ids":[],"tissue_ontology_term_ids":["UBERON:0002113"]},"is_rollup":true}'
# 
# h <- new_handle(url = url)
# handle_setheaders(h, headers)
# handle_setopt(h, customrequest = "POST")
# handle_setopt(h, postfields = body)
# handle_setopt(h, verbose = TRUE)
# 
# result <- tryCatch({
#   curl_fetch_memory(handle = h, url=url)
# }, error = function(e) {
#   stop(paste("An error occurred:", conditionMessage(e)))
# })
# 
# if (!inherits(result, "try-error")) {
#   response <- parse_memory(result)
#   content <- rawToChar(response$content)
#   print(content)
# }
# 
# handle_close(h)
# 
# 
# 
# #####
# library(httr)
# library(jsonlite)
# 
# 
# 
# # Create the body of the POST request
# json_body <- list(
#   "filter" = list(
#     "dataset_ids" = list("0b4a15a7-4e9e-4555-9733-2423e5c66469"),
#     "gene_ontology_term_ids" = list("ENSG00000089597"))
# )
# 
# 
# request <- POST(
#   url = input_url, 
#   # add_headers("Authorization" = paste("Bearer", notion_api_key),
#   #             "Notion-Version" = "2021-05-13"),
#   body = json_body, 
#   encode = "json"
# )
# 
# 
# 
# 
# # curl -v -H "Content-Type: application/json" -X POST -d  
# # '{"filter":{"dataset_ids":["0b4a15a7-4e9e-4555-9733-2423e5c66469"],"development_stage_ontology_term_ids":[],"disease_ontology_term_ids":[],"gene_ontology_term_ids":["ENSG00000089597"],
# # "organism_ontology_term_id":"NCBITaxon:9606","self_reported_ethnicity_ontology_term_ids":[],"sex_ontology_term_ids":[],"tissue_ontology_term_ids":["UBERON:0002113"]},"is_rollup":true}' 
# # https://api.cellxgene.cziscience.com/wmg/v1/query
# # 
# 
# curl -X POST 'https://api.notion.com/v1/databases/897e5a76ae524b489fdfe71f5945d1af' \
# -H 'Authorization: Bearer '"$NOTION_API_KEY"'' \
# -H 'Notion-Version: 2021-05-13' \
# -H "Content-Type: application/json" \
# --data '{
#       "filter": {
#         "or": [
#           {
#             "property": "In stock",
#                     "checkbox": {
#                         "equals": true
#                     }
#           },
#           {
#                     "property": "Cost of next trip",
#                     "number": {
#                         "greater_than_or_equal_to": 2
#                     }
#                 }
#             ]
#         },
#       "sorts": [
#         {
#           "property": "Last ordered",
#           "direction": "ascending"
#         }
#       ]
#     }'
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Create the body of the POST request
# json_body <- list(
#   "filter" = list(
#     "or" = list(
#       list(
#         "property" = "In stock",
#         "checkbox" = list(
#           "equals" = "true"
#         )
#       ),
#       list(
#         "property" = "Cost of next trip",
#         "number" = list(
#           "greater_than_or_equal_to" = 2
#         )
#       )
#     )
#   ),
#   "sorts" = list(
#     list(
#       "property" = "Last ordered",
#       "direction" = "ascending"
#     )
#   )
# )
# 
# # Make post request
# request <- POST(
#   url = "https://api.notion.com/v1/databases/897e5a76ae524b489fdfe71f5945d1af", 
#   add_headers("Authorization" = paste("Bearer", notion_api_key),
#               "Notion-Version" = "2021-05-13"),
#   body = json_body, 
#   encode = "json"
# )
# 
# 
#     
#     
#     








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
# x <- read_json(content)

# https://www.geeksforgeeks.org/how-to-read-large-json-file-in-r/

setwd("/Users/nrank/Desktop/BioInf/halbritter/nephro_gene_score/")

writeLines(content[[1]], "junk_content.txt")                 # Apply writeLines function

data <- jsonlite::read_json("junk_content.txt")

head(data, 3)



