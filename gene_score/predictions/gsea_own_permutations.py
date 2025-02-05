# %load_ext autoreload
# %autoreload 2
# import basic modules
from datetime import datetime
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
import requests
import sys
import urllib.request
import yaml



# import third-party modules
import gseapy as gp

# set options
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"
script_path = "/gene_score/"

# read configs
CONFIG_FILE = os.getenv('CONFIG_FILE')

with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}{script_path}")

subset = 4  # CAVE: change np.random.seed with every subset
print(f"subset: {subset}")

# TODO: customize results version 
ID = 97
results_date = "2024-03-22"

###############################################################################################
# get GO terms and associated genes from library 'GO_Biological_Process_2023'
GO_BP_human_lib = gp.get_library("GO_Biological_Process_2023", organism='Human')

# create a df of GO terms and associated genes
GO_BP_human = pd.DataFrame(list(GO_BP_human_lib.items()), columns=['full_GO_term', 'genes'])

# extract GO_term
GO_BP_human['GO_term'] = GO_BP_human['full_GO_term'].str.extract(r'\((GO:\d+)\)')
###############################################################################################


###############################################################################################
### check if all genes of 'GO_Biological_Process_2023' have the same notation as in NGS

# get all genes of 'GO_Biological_Process_2023'
GO_BP_human_exploded = GO_BP_human['genes'].explode()
GO_BP_human_all_genes = list(set(GO_BP_human_exploded.tolist()))

## check which genes have no match in NGS
# load Nephro Gene Score (NGS)
NGS = pd.read_csv(f"predictions/results/NGS_predictions_ID{ID}_all_{results_date}.csv.gz")
NGS['symbol'] = NGS['symbol'].str.upper()
NGS['hgnc_id'] = "HGNC:" + NGS['hgnc_id_int'].astype(str)

unmatched_genes = [i for i in GO_BP_human_all_genes if i not in NGS['symbol'].values]

## create a mapping form unmatched genes to HGNC IDs

# download HGNC annotated table from kidney-genetics
hgnc_annotated_url = f"https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/B_AnnotationHGNC/results/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version']}.csv.gz"
hgnc_annotated_dest_file = f"raw/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version']}.csv.gz"

# check if the file already exists
if not os.path.exists(hgnc_annotated_dest_file):
    # download the file
    urllib.request.urlretrieve(hgnc_annotated_url, hgnc_annotated_dest_file)
    print(f"The file '{hgnc_annotated_dest_file}' has been downloaded.")
else:
    print(f"The file '{hgnc_annotated_dest_file}' already exists. Skipping the download.")

HGNC_table = pd.read_csv(f"raw/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version']}.csv.gz", low_memory=False)


def get_hgnc_id_from_alias_or_prev_symbol(symbol, HGNC_table=None):
    """
    Function that returns HGNC ID from symbol (via alias or previous symbol).
    Example: get_hgnc_id_from_alias_or_prev_symbol('PKD1')
    """
    
    if HGNC_table is None:
        HGNC_table = pd.read_csv(f"raw/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version']}.csv.gz", low_memory=False)
        
    HGNC_sub = HGNC_table[['hgnc_id', 'alias_symbol', 'prev_symbol']].copy()

    # fill NaN values with empty strings
    HGNC_sub['alias_symbol'] = HGNC_sub['alias_symbol'].fillna('')
    HGNC_sub['prev_symbol'] = HGNC_sub['prev_symbol'].fillna('')

    # combine alias symbols and previous symbols
    HGNC_sub['alias_or_previous'] = HGNC_sub['alias_symbol'] + "|" + HGNC_sub['prev_symbol'] 
    HGNC_sub['alias_or_previous'] = HGNC_sub['alias_or_previous'].str.upper()
    HGNC_sub = HGNC_sub.query("alias_or_previous != '|'")

    # split values in 'alias_symbol' column by '|', then explode into separate rows
    HGNC_sub_exploded = HGNC_sub.assign(alias_or_previous=HGNC_sub['alias_or_previous'].str.split('|')).explode('alias_or_previous')
    
    # filter rows where 'alias_or_previous' contains symbol
    symbol_rows = HGNC_sub_exploded[HGNC_sub_exploded['alias_or_previous'] == symbol]
    
    # check if multiple HGNC IDs have the same alias/previous symbol
    if len(symbol_rows) > 1:
        raise ValueError(f"Multiple HGNC IDs for {symbol}.")
    
    # if no alias/previous symbol available, return None
    if len(symbol_rows) == 0:
        return None
        
    # get HGNC ID
    hgnc_id_symbol = symbol_rows['hgnc_id'].iloc[0]
    
    return hgnc_id_symbol


gene_mapping = {}

for gene in unmatched_genes:
    mapped_hgnc_id = get_hgnc_id_from_alias_or_prev_symbol(symbol=gene, HGNC_table=HGNC_table)
    
    if mapped_hgnc_id is not None:
        gene_mapping[mapped_hgnc_id] = gene
    

# function that returns the value of gene_mapping if hgnc_id is in gene_mapping, else the original symbol
def convert_symbol(row):
    if row['hgnc_id'] in gene_mapping:
        return gene_mapping[row['hgnc_id']]
    else:
        return row['symbol']

# apply the function to create the new column 'symbol_converted'
NGS['symbol_converted'] = NGS.apply(lambda row: convert_symbol(row), axis=1)
###############################################################################################


###############################################################################################
## get kidney GO terms

def get_descendants(GO_term):
    """Function that returns all descendants of a GO term using the QuickGO API"""
    
    # extract number of GO term
    GO_number = GO_term.split(":")[-1]

    # set request URL
    requestURL = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A{GO_number}/descendants?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"

    # get result
    res = requests.get(requestURL, headers={ "Accept" : "application/json"})
    
    # parse text
    responseBody = json.loads(res.text)
    
    # access descendants
    descendants = responseBody['results'][0]['descendants']

    return descendants


# get descendants of GO:0072001 ("renal system development")
desc_GO_0072001 = get_descendants("GO:0072001")

# get descendants of GO:0003014 ("renal system process")
desc_GO_0003014 = get_descendants("GO:0003014")

# combine all kidney GO terms
kidney_GO_terms = list(set(["GO:0072001"] + desc_GO_0072001 + ["GO:0003014"] + desc_GO_0003014))    
    
###############################################################################################
# scale NGS
scale = False # if True, NGS will be scaled from -1 to +1

# rank df by NGS value
NGS_ranked = NGS.sort_values('NGS', ascending = False).reset_index(drop = True)[['symbol_converted', 'NGS']]

# scale NGS value so that it reaches from -1 to +1
if scale:
#     NGS_ranked['NGS'] = 2 * NGS_ranked['NGS'] - 1
    NGS_ranked['NGS'] = NGS_ranked['NGS'] - 0.5
    
    
###############################################################################################
    
###############################################################################################    
## GSEA analysis with ranked gene list (ranked by NGS)    
# perform GSEA
pre_res = gp.prerank(rnk = NGS_ranked,
                     gene_sets = ['GO_Biological_Process_2023'], 
                     seed = 42, 
                     permutation_num = 1000,
                     outdir = None,
                     min_size = 0,
                     max_size = 5000000 # all gene sets included for max_size = 5000000 and min_size = 0
                    )

# get results df
NGS_pre_res = pre_res.res2d

# create a new column with GO term
NGS_pre_res['GO_term'] = NGS_pre_res['Term'].str.extract(r'\((GO:\d+)\)')

# create a new column that indicates if GO term is in kidney GO terms
NGS_pre_res['is_kidney_GO'] = NGS_pre_res['GO_term'].isin(kidney_GO_terms)
###############################################################################################


###############################################################################################
## Enrichment analysis for kidney GO terms in GO terms ranked by NES of previous analysis
# Method: as in GSEA with GO_terms as "genes" and kidney GO terms as "gene set"

# sort results of previous GSEA by Normalized enrichment score (NES)
NGS_pre_res_NES_ranked = NGS_pre_res.sort_values(by='NES', ascending=False).reset_index()

NGS_pre_res_NES_ranked['NES'] = NGS_pre_res_NES_ranked['NES'].astype(float)

# create dictionary of kidney GO terms as input for gp.prerank()
set_kidney_GO = {'kid_GO_terms':[i for i in kidney_GO_terms if i in NGS_pre_res['GO_term'].values]}

# perform analysis
pre_res_kid_GO = gp.prerank(rnk = NGS_pre_res_NES_ranked[['GO_term', 'NES']],
                            gene_sets = set_kidney_GO,
                            background = None,
                            seed = 42, 
                            permutation_num = 1000,
                            outdir = None,
                            min_size = 0,
                            max_size = 50000 # corresponds to no limitation on set size
                           )

# get results
# pre_res_kid_GO.res2d

current_date = datetime.now().strftime("%Y-%m-%d")

pre_res_kid_GO.res2d.to_csv(f"ID{ID}_pre_res_kid_GO_real_{current_date}.csv", index=False)


set_kidney_GO = {'kid_GO_terms':[i for i in kidney_GO_terms if i in NGS_pre_res['GO_term'].values]}



def gsea_manual_perm(rnk, set_kidney_GO):
    ## GSEA analysis with ranked gene list (ranked by NGS)    
    # perform GSEA
    
    pre_res = gp.prerank(rnk = rnk,
                     gene_sets = ['GO_Biological_Process_2023'], 
                     seed = 42, 
                     permutation_num = 1000,
                     outdir = None,
                     min_size = 0,
                     max_size = 5000000 # all gene sets included for max_size = 5000000 and min_size = 0
                    )

    # get results df
    NGS_pre_res = pre_res.res2d

    # create a new column with GO term
    NGS_pre_res['GO_term'] = NGS_pre_res['Term'].str.extract(r'\((GO:\d+)\)')

#     # create a new column that indicates if GO term is in kidney GO terms
#     NGS_pre_res['is_kidney_GO'] = NGS_pre_res['GO_term'].isin(kidney_GO_terms)
    
    ## Enrichment analysis for kidney GO terms in GO terms ranked by NES of previous analysis
    # Method: as in GSEA with GO_terms as "genes" and kidney GO terms as "gene set"

    # sort results of previous GSEA by Normalized enrichment score (NES)
    NGS_pre_res_NES_ranked = NGS_pre_res.sort_values(by='NES', ascending=False).reset_index()

    NGS_pre_res_NES_ranked['NES'] = NGS_pre_res_NES_ranked['NES'].astype(float)

    # create dictionary of kidney GO terms as input for gp.prerank()
#     set_kidney_GO = {'kid_GO_terms':[i for i in kidney_GO_terms if i in NGS_pre_res['GO_term'].values]}

    # perform analysis
    pre_res_kid_GO = gp.prerank(rnk = NGS_pre_res_NES_ranked[['GO_term', 'NES']],
                                gene_sets = set_kidney_GO,
                                background = None,
                                seed = 42, 
                                permutation_num = 1000,
                                outdir = None,
                                min_size = 0,
                                max_size = 50000 # corresponds to no limitation on set size
                               )
    
    # get results
    return pre_res_kid_GO.res2d


pre_res_kid_GO_permutations = pd.DataFrame()


for i in np.arange(250) + 750:
    print(i)
    sys.stdout.flush() 
    np.random.seed(i)  # CAVE!! change seed with other subsets!!
    
    # shuffle NGS ranking
    NGS_shuffled = NGS_ranked.copy()
    NGS_shuffled['symbol_converted'] = np.random.permutation(NGS_shuffled['symbol_converted'].values)
    
    new_row = gsea_manual_perm(rnk=NGS_shuffled, set_kidney_GO=set_kidney_GO)
    
    pre_res_kid_GO_permutations = pd.concat([pre_res_kid_GO_permutations, new_row])

        
    
pre_res_kid_GO_permutations.to_csv(f"ID{ID}_pre_res_kid_GO_permutations_subset{subset}_raw_NGS_{current_date}.csv", index=False)




