# Title: helper functions

# function for left join
# create a left_join function that joins by second identifier, if first identifier does not hit
left_join_rescue <- function(x, y, by1, by2) {
  
  # try to join by first column ('by1')
  y_sub1 <- y %>% dplyr::select(-{{ by2 }})
  by_by1_match <- inner_join(x, y_sub1, by = by1) # joined df that was joined by column 'by1'
  by_by1_no_match <- anti_join(x, y_sub1, by = by1) # df that had no match by column 'by1'
  
  # try to join by first column ('by2')
  y_sub2 <- y %>% dplyr::select(-{{ by1 }})
  by_by2_match <- inner_join(by_by1_no_match, y_sub2, by = by2) # joined df that was joined by column 'by1'
  no_match <- anti_join(by_by1_no_match, y_sub2, by = by2) # df that had no match at all (neither by 'by1' nor by 'by2')
  
  # combine subset dataframes
  res <- dplyr::bind_rows(by_by1_match, by_by2_match, no_match)
  
  return(res)
}

# create a left_join function that joins by hgnc_id (from symbol), if first identifier does not hit
left_join_rescue_symbol <- function(x, y, by1) {
  
  # try to join by first column ('by1')
  y_sub1 <- y %>% dplyr::select(-symbol)
  by_by1_match <- inner_join(x, y_sub1, by = by1) # joined df that was joined by column 'by1'
  by_by1_no_match_x <- anti_join(x, y_sub1, by = by1) # subdf of x that had no match by column 'by1'
  
  by_by1_no_match_y <- anti_join(y, x, by = by1) %>% filter(symbol %in% all_prot_coding_gene_symbols)
  
  # get hgnc_id for symbols
  by_by1_no_match_y$hgnc_id_int <- hgnc_id_from_symbol_grouped(tibble(value = by_by1_no_match_y$symbol)) 
  
  by_by1_no_match_y <- by_by1_no_match_y %>% dplyr::select(-{{ by1 }}, -symbol)
  
  # in case two symbols match the same hgnc_id use only the first row
  by_by1_no_match_y <- by_by1_no_match_y[!duplicated(by_by1_no_match_y$hgnc_id_int), ]
  
  # try to join by "hgnc_id_int"
  by_hgnc_id_match <- inner_join(by_by1_no_match_x, by_by1_no_match_y, by = "hgnc_id_int") 
  no_match <- anti_join(by_by1_no_match_x, by_by1_no_match_y, by = "hgnc_id_int") # df that had no match at all 
  
  # combine subset dataframes
  res <- dplyr::bind_rows(by_by1_match, by_hgnc_id_match, no_match)
  
  return(res)
}

# NOTE: these functions work slightly different, in both cases matching issues arrive. In general left_join_rescue_symbol() should be used.

