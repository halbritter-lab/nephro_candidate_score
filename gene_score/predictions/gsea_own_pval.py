#!/usr/bin/env python
# coding: utf-8



### GSEA analysis #### - converted from .ipynb

# %load_ext autoreload
# %autoreload 2
# import basic modules
from datetime import datetime
import json
import matplotlib.pyplot as plt
import pandas as pd
import random
import requests
import sys
import urllib.request
import yaml
import pickle
import os

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


# In[4]:


# Notes: 
# https://www.youtube.com/watch?v=Yi4d7JIlAsM
# https://www.cs.tufts.edu/comp/167/gsea.pdf for Gene % explanation
# supplement of https://www.pnas.org/doi/10.1073/pnas.0506580102 => original GSEA method, calculation of ES


# In[5]:


# Notes on gp.prerank():
# - increase of permutation_num => changes NES (+ order of GO terms), FDR-q, FWER, but NOT ES, Tag, Gene 
# - limitation of max_size => analysis is only done for GO-terms with <= max_size genes 
# => ES, NES, order stays the same, but FDR-q, FWER change as fewer pathways are analysed (fewer multiple testing)
# - alternation of NGS changes ES etc. multiplication doesn't change it, but addition of certain value changes it
# => rank alone is not criteria, but relation between the ranks??

# - max size of gene sets in GO_BP is 2002 genes # set([int(i.split("/")[-1]) for i in pre_res.res2d['Tag %'].values])


# In[6]:


# show gseapy libraries for organism = 'Human'
gp.get_library_name(organism = "Human")


# In[7]:


# get GO terms and associated genes from library 'GO_Biological_Process_2023'
GO_BP_human_lib = gp.get_library("GO_Biological_Process_2023", organism='Human')

# create a df of GO terms and associated genes
GO_BP_human = pd.DataFrame(list(GO_BP_human_lib.items()), columns=['full_GO_term', 'genes'])

# extract GO_term
GO_BP_human['GO_term'] = GO_BP_human['full_GO_term'].str.extract(r'\((GO:\d+)\)')

# # get GO terms and associated genes from library 'GO_Cellular_Component_2023'
# GO_CC_human_lib = gp.get_library("GO_Cellular_Component_2023", organism='Human')
# GO_CC_human = pd.DataFrame(list(GO_CC_human_lib.items()), columns=['full_GO_term', 'genes'])
# GO_CC_human['GO_term'] = GO_CC_human['full_GO_term'].str.extract(r'\((GO:\d+)\)')

# Note: no kidney GO terms in 'GO_Cellular_Component_2023'


# # get GO terms and associated genes from library 'GO_Molecular_Function_2023'
# GO_MF_human_lib = gp.get_library("GO_Molecular_Function_2023", organism='Human')
# GO_MF_human = pd.DataFrame(list(GO_MF_human_lib.items()), columns=['full_GO_term', 'genes'])
# GO_MF_human['GO_term'] = GO_MF_human['full_GO_term'].str.extract(r'\((GO:\d+)\)')

# Note: no kidney GO terms in 'GO_Molecular_Function_2023'


# In[8]:


### check if all genes of 'GO_Biological_Process_2023' have the same notation as in NGS

# get all genes of 'GO_Biological_Process_2023'
GO_BP_human_exploded = GO_BP_human['genes'].explode()
GO_BP_human_all_genes = list(set(GO_BP_human_exploded.tolist()))

## check which genes have no match in NGS
# load Nephro Gene Score (NGS)
NGS = pd.read_csv(f"predictions/results/NGS_predictions_ID{config_vars['ID_ngs']}_all_{config_vars['ngs_results_date']}.csv.gz")
NGS['symbol'] = NGS['symbol'].str.upper()
NGS['hgnc_id'] = "HGNC:" + NGS['hgnc_id_int'].astype(str)

unmatched_genes = [i for i in GO_BP_human_all_genes if i not in NGS['symbol'].values]

## create a mapping form unmatched genes to HGNC IDs

# download HGNC annotated table from kidney-genetics
hgnc_annotated_url = f"https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/B_AnnotationHGNC/results/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version_gs']}.csv.gz"
hgnc_annotated_dest_file = f"raw/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version_gs']}.csv.gz"

# check if the file already exists
if not os.path.exists(hgnc_annotated_dest_file):
    # download the file
    urllib.request.urlretrieve(hgnc_annotated_url, hgnc_annotated_dest_file)
    print(f"The file '{hgnc_annotated_dest_file}' has been downloaded.")
else:
    print(f"The file '{hgnc_annotated_dest_file}' already exists. Skipping the download.")

HGNC_table = pd.read_csv(f"raw/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version_gs']}.csv.gz", low_memory=False)


def get_hgnc_id_from_alias_or_prev_symbol(symbol, HGNC_table=None):
    """
    Function that returns HGNC ID from symbol (via alias or previous symbol).
    Example: get_hgnc_id_from_alias_or_prev_symbol('PKD1')
    """
    
    if HGNC_table is None:
        HGNC_table = pd.read_csv(f"raw/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version_gs']}.csv.gz", low_memory=False)
        
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



# In[ ]:





# In[9]:


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
    
# save as csv
current_date = datetime.now().strftime("%Y-%m-%d")
# pd.DataFrame({'kidney_GO_term': kidney_GO_terms}).to_csv(f"predictions/results/kidney_GO_terms_{current_date}.csv.gz", index=False, compression='gzip')


# In[12]:


# get number of genes that are present both in GO_BP and kidney_GO_terms
common_terms = set(GO_BP_human['GO_term']).intersection(kidney_GO_terms)
len(common_terms)


# In[13]:


## GSEA analysis with ranked gene list (ranked by NGS)
scale = False # if True, NGS will be scaled from -1 to +1

# rank df by NGS value
NGS_ranked = NGS.sort_values('NGS', ascending = False).reset_index(drop = True)[['symbol_converted', 'NGS']]

# scale NGS value so that it reaches from -1 to +1
if scale:
    NGS_ranked['NGS'] = 2 * NGS_ranked['NGS'] - 1

# perform GSEA
pre_res = gp.prerank(rnk = NGS_ranked,
                     gene_sets = ['GO_Biological_Process_2023'], 
                     seed = 42, 
                     permutation_num = 1000, # CHANGE!!!!
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

# # save results
current_date = datetime.now().strftime("%Y-%m-%d")


if scale:
    NGS_pre_res.to_csv(f"predictions/results/GSEA1_NGS_scaled_ranking_{current_date}.csv.gz", index=False, compression='gzip')
    with open(f"predictions/results/GSEA1_pre_res_scaled_ranking_{current_date}.pkl", 'wb') as f:
        pickle.dump(pre_res, f)
else:
    NGS_pre_res.to_csv(f"predictions/results/GSEA1_NGS_ranking_{current_date}.csv.gz", index=False, compression='gzip')
    with open(f"predictions/results/GSEA1_pre_res_ranking_{current_date}.pkl", 'wb') as f:
        pickle.dump(pre_res, f)


# In[25]:


# save only top 20
gsea1_top20 = pd.read_csv(f"predictions/results/GSEA_NGS_ranking_{current_date}.csv.gz")[['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val', 'FWER p-val',
       'Tag %', 'is_kidney_GO']].head(20)

gsea1_top20['Process'] = [i.split("__")[1] for i in gsea1_top20['Term'].values]

gsea1_top20 = gsea1_top20[['Process', 'ES', 'NES', 'NOM p-val', 'FDR q-val', 'FWER p-val','Tag %', 'is_kidney_GO']]
gsea1_top20[['ES', 'NES', 'NOM p-val', 'FDR q-val', 'FWER p-val']] = gsea1_top20[['ES', 'NES', 'NOM p-val', 'FDR q-val', 'FWER p-val']].round(3)

gsea1_top20.to_csv(f"predictions/results/GSEA_NGS_ranking_top20_{current_date}.csv", index=False)


# In[12]:


# Explanations of columns

# ES (Enrichment Score): The Enrichment Score reflects the degree to which a gene set is overrepresented at the extremes (top or bottom) of the ranked list of genes. It measures the enrichment of genes associated with a particular set of functions (e.g., kidney GO terms) towards one end of the ranked list. Higher positive or negative ES values indicate stronger enrichment towards the top or bottom of the ranked list, respectively.

# NES (Normalized Enrichment Score): The Normalized Enrichment Score is the primary statistic used to rank and compare gene sets in GSEA. It's the ES normalized to account for differences in gene set size and correlation structure. NES values provide a more accurate assessment of the enrichment significance compared to the ES alone.

# NOM p-val (Nominal p-value): The Nominal p-value represents the statistical significance of the enrichment score. It indicates the probability of observing the given NES by chance, without correction for multiple hypothesis testing.

# FDR q-val (False Discovery Rate q-value): The False Discovery Rate (FDR) q-value is a measure that adjusts the nominal p-values for multiple hypothesis testing using methods such as the Benjamini-Hochberg procedure. It represents the expected proportion of false positives among the tests that are called significant. Lower FDR q-values indicate stronger evidence for the significance of the enrichment.

# FWER p-val (Family-Wise Error Rate p-value): The Family-Wise Error Rate (FWER) p-value is another measure of statistical significance that controls for multiple testing. It represents the probability of obtaining the observed NES or a more extreme value under the assumption of no association between the gene set and the phenotype. FWER p-values are typically more conservative than FDR q-values.

# Tag % (Tag percentage): The Tag percentage represents the percentage of genes in the input list (ranked list of all human genes) that overlap with the genes in the gene set of interest (e.g., kidney GO terms). It indicates the extent to which the gene set is represented in your input data.

# Gene % (Gene percentage): The Gene percentage represents the percentage of genes in the gene set of interest (e.g., kidney GO terms) that are found in the input list. It indicates how well the genes from your gene set are represented in your input data.


# In[25]:


# Plots
if scale:
    with open(f"predictions/results/GSEA1_pre_res_scaled_ranking_{current_date}.pkl", 'rb') as f:
        pre_res = pickle.load(f)
else: 
    with open(f"predictions/results/GSEA1_pre_res_ranking_{current_date}.pkl", 'rb') as f:
        pre_res = pickle.load(f)

terms = pre_res.res2d.Term
axs = pre_res.plot(terms=terms[0:5],
                   #legend_kws={'loc': (1.2, 0)}, # set the legend loc
                   show_ranking=True, # whether to show the second yaxis
                   figsize=(3,4)
                  )
# or use this to have more control on the plot
# from gseapy import gseaplot2
# terms = pre_res.res2d.Term[1:5]
# hits = [pre_res.results[t]['hits'] for t in terms]
# runes = [pre_res.results[t]['RES'] for t in terms]
# fig = gseaplot2(terms=terms, ress=runes, hits=hits,
#               rank_metric=gs_res.ranking,
#               legend_kws={'loc': (1.2, 0)}, # set the legend loc
#               figsize=(4,5)) # rank_metric=pre_res.ranking


# In[17]:


# ## Enrichment analysis for kidney GO terms in GO terms ranked by NES of previous analysis
# # Method: as in GSEA with GO_terms as "genes" and kidney GO terms as "gene set"

# # sort results of previous GSEA by Normalized enrichment score (NES)
# NGS_pre_res_NES_ranked = NGS_pre_res.sort_values(by='NES', ascending=False).reset_index()

# NGS_pre_res_NES_ranked['NES'] = NGS_pre_res_NES_ranked['NES'].astype(float)

# # create dictionary of kidney GO terms as input for gp.prerank()
# set_kidney_GO = {'kid_GO_terms':[i for i in kidney_GO_terms if i in NGS_pre_res['GO_term'].values]}

# # perform analysis
# pre_res_kid_GO = gp.prerank(rnk = NGS_pre_res_NES_ranked[['GO_term', 'NES']],
#                             gene_sets = set_kidney_GO,
#                             background = None,
#                             seed = 42, 
#                             permutation_num = 1000,
#                             outdir = None,
#                             min_size = 0,
#                             max_size = 50000 # corresponds to no limitation on set size
#                            )

# # get results
# pre_res_kid_GO.res2d


# In[18]:


# # Plot result 
# terms = pre_res_kid_GO.res2d.Term

# axs = pre_res_kid_GO.plot(terms=terms[0],
#                    #legend_kws={'loc': (1.2, 0)}, # set the legend loc
#                    show_ranking=True, 
#                    figsize=(15,10)
#                   )


# In[17]:


# ## Enrichment analysis for kidney GO terms in GO terms ranked by NES of previous analysis

# # Problem: only checks if there are more elements of the unrankedSublist in the upper
# # n = len(unrankedSublist) than by chance. Exact ranking not relevant.

# def test_statistic(ranked_list, unranked_sublist, occur_list):
#     # consider top portion of ranked list
#     top_portion = ranked_list[:len(unranked_sublist)]  
    
#     # count occurences of elements of unrankedSublist in the top portion of the ranked list
#     occurrences = sum(1 for element in unranked_sublist if element in top_portion)
    
#     # append occurences list
#     occur_list.append(occurrences)
#     return occurrences, occur_list

# def permute(ranked_list):
#     # create a copy of the ranked list
#     permuted_ranked_list = ranked_list[:] 
    
#     # shuffle the copied list
#     random.shuffle(permuted_ranked_list)  
#     return permuted_ranked_list

# def permutation_test(ranked_list, unranked_sublist, num_permutations, occur_list):
#     random.seed(1)
#     observed_statistic, occur_list = test_statistic(ranked_list, unranked_sublist, occur_list)
#     permuted_statistics = []
    
#     for _ in range(num_permutations):
#         permuted_ranked_list = permute(ranked_list)  # Randomly permute ranked list
#         permuted_statistic, occur_list = test_statistic(permuted_ranked_list, unranked_sublist, occur_list)
#         permuted_statistics.append(permuted_statistic)
    
#     # calculate p-value
#     p_value = sum(1 for stat in permuted_statistics if stat >= observed_statistic) / num_permutations
#     return p_value, occur_list



# ranked_list = NGS_pre_res.sort_values(by='NES', ascending=False)['GO_term'].copy().values
# unranked_sublist = [i for i in kidney_GO_terms if i in ranked_list]


# p_value, occur_list = permutation_test(ranked_list = ranked_list,
#                                        unranked_sublist = unranked_sublist,
#                                        num_permutations = 1000, 
#                                        occur_list=[]
#                                       )

# print("p-value:", p_value)


# In[ ]:


# Finally chosen method: shuffle NGS randomly, repeat GSEA1, and GSEA2, do this 1000 times and get p-value by
# comparing the NES distribution to the real observed NES
# see scripts: gsea_own_permutations.py, gsea_own_permutations.sh

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


gsea_manual_perm(rnk = NGS_ranked, set_kidney_GO=set_kidney_GO)


# In[10]:


get_ipython().system('pwd')


# In[13]:


# compare distribution of GSEA2 ES values yielded by random shuffling of N-GS in GSEA1 with the observed ("real") value

subset1 = pd.read_csv("predictions/results/ID97_pre_res_kid_GO_permutations_subset1_raw_NGS_2024-06-03.csv.gz")
subset2 = pd.read_csv("predictions/results/ID97_pre_res_kid_GO_permutations_subset2_raw_NGS_2024-06-03.csv.gz")
subset3 = pd.read_csv("predictions/results/ID97_pre_res_kid_GO_permutations_subset3_raw_NGS_2024-06-03.csv.gz")
subset4 = pd.read_csv("predictions/results/ID97_pre_res_kid_GO_permutations_subset4_raw_NGS_2024-06-03.csv.gz")

all_perm = pd.concat([subset1, subset2, subset3, subset4])
all_perm
plt.hist(all_perm['ES'], bins=100)

real_NGS_res =  pd.read_csv("predictions/results/ID97_pre_res_kid_GO_real_2024-06-03.csv")
real_NGS_res


# In[14]:


import numpy as np
np.sum([1 for i in all_perm.ES.values if i >= real_NGS_res.NES.values])/1000


# In[15]:


# compare distribution of GSEA2 NES values yielded by random shuffling of N-GS in GSEA1 with the observed ("real") value
subset1 = pd.read_csv("predictions/results/ID97_pre_res_kid_GO_permutations_subset1_raw_NGS_2024-06-03.csv.gz")
subset2 = pd.read_csv("predictions/results/ID97_pre_res_kid_GO_permutations_subset2_raw_NGS_2024-06-03.csv.gz")
subset3 = pd.read_csv("predictions/results/ID97_pre_res_kid_GO_permutations_subset3_raw_NGS_2024-06-03.csv.gz")
subset4 = pd.read_csv("predictions/results/ID97_pre_res_kid_GO_permutations_subset4_raw_NGS_2024-06-03.csv.gz")

all_perm = pd.concat([subset1, subset2, subset3, subset4])
all_perm
plt.hist(all_perm['NES'], bins=100)

real_NGS_res =  pd.read_csv("predictions/results/ID97_pre_res_kid_GO_real_2024-06-03.csv")
real_NGS_res


# In[16]:


import numpy as np
np.sum([1 for i in all_perm.NES.values if i >= real_NGS_res.NES.values])/1000


# In[ ]:




