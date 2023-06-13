# Title: helper functions for getting features 

# load libraries
library(tidyverse)
library(jsonlite)

# functions to retrieve children from an HPO term
HPO_children_from_term <- function(term_input_id) {
  hpo_term_response <- fromJSON(paste0("https://hpo.jax.org/api/hpo/term/", URLencode(term_input_id, reserved=T)))
  hpo_term_children <- as_tibble(hpo_term_response$relations$children)
  return(hpo_term_children)
}

HPO_all_children_from_term <- function(term_input) {
  children_list <- HPO_children_from_term(term_input)
  all_children_list <<- append(all_children_list, term_input)
  if(length(children_list)!=0)
  {
    for (p in children_list$ontologyId) {
      all_children_list <<- append(all_children_list, p)
      Recall(p)
    }
  }
  all_children_tibble <- as_tibble(unlist(all_children_list)) %>% unique
  return(all_children_tibble)
}

