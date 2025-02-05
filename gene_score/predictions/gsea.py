#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.environ['CONFIG_FILE'] = '/fast/work/users/rankn_c/halbritter/nephro_candidate_score/gene_score/training/config_NCS.yml' # TODO: change


# In[9]:


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


# In[3]:


# Notes: 
# https://www.youtube.com/watch?v=Yi4d7JIlAsM
# https://www.cs.tufts.edu/comp/167/gsea.pdf for Gene % explanation
# supplement of https://www.pnas.org/doi/10.1073/pnas.0506580102 => original GSEA method, calculation of ES


# In[4]:


# TODO: customize results version 
ID = 97
results_date = "2024-03-22"


# CONTINUE HERE: 
# - check if NGS symbol notation is the same as in enrichr libraries

# Notes on gp.prerank():
# - increase of permutation_num => changes NES (+ order of GO terms), FDR-q, FWER, but NOT ES, Tag, Gene 
# - limitation of max_size => analysis is only done for GO-terms with <= max_size genes 
# => ES, NES, order stays the same, but FDR-q, FWER change as fewer pathways are analysed (fewer multiple testing)
# - alternation of NGS changes ES etc. multiplication doesn't change it, but addition of certain value changes it
# => rank alone is not criteria, but relation between the ranks??

# - max size of gene sets in GO_BP is 2002 genes # set([int(i.split("/")[-1]) for i in pre_res.res2d['Tag %'].values])


# In[5]:


# show gseapy libraries for organism = 'Human'
gp.get_library_name(organism = "Human")


# In[ ]:





# In[6]:


# get GO terms and associated genes from library 'GO_Biological_Process_2023'
GO_BP_human_lib = gp.get_library("GO_Biological_Process_2023", organism='Human')

# create a df of GO terms and associated genes
GO_BP_human = pd.DataFrame(list(GO_BP_human_lib.items()), columns=['full_GO_term', 'genes'])

# extract GO_term
GO_BP_human['GO_term'] = GO_BP_human['full_GO_term'].str.extract(r'\((GO:\d+)\)')

GO_BP_human.head(5)


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


# In[7]:


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



# In[10]:


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
pd.DataFrame({'kidney_GO_term': kidney_GO_terms}).to_csv(f"predictions/results/kidney_GO_terms_{current_date}.csv.gz", index=False, compression='gzip')


# In[11]:


## GSEA analysis with ranked gene list (ranked by NGS)

# rank df by NGS value
NGS_ranked = NGS.sort_values('NGS', ascending = False).reset_index(drop = True)[['symbol_converted', 'NGS']]

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

# save results
NGS_pre_res.to_csv(f"predictions/results/GSEA_NGS_ranking_{current_date}.csv.gz", index=False, compression='gzip')


# In[12]:


# Explanations of columns

# ES (Enrichment Score): The Enrichment Score reflects the degree to which a gene set is overrepresented at the extremes (top or bottom) of the ranked list of genes. It measures the enrichment of genes associated with a particular set of functions (e.g., kidney GO terms) towards one end of the ranked list. Higher positive or negative ES values indicate stronger enrichment towards the top or bottom of the ranked list, respectively.

# NES (Normalized Enrichment Score): The Normalized Enrichment Score is the primary statistic used to rank and compare gene sets in GSEA. It's the ES normalized to account for differences in gene set size and correlation structure. NES values provide a more accurate assessment of the enrichment significance compared to the ES alone.

# NOM p-val (Nominal p-value): The Nominal p-value represents the statistical significance of the enrichment score. It indicates the probability of observing the given NES by chance, without correction for multiple hypothesis testing.

# FDR q-val (False Discovery Rate q-value): The False Discovery Rate (FDR) q-value is a measure that adjusts the nominal p-values for multiple hypothesis testing using methods such as the Benjamini-Hochberg procedure. It represents the expected proportion of false positives among the tests that are called significant. Lower FDR q-values indicate stronger evidence for the significance of the enrichment.

# FWER p-val (Family-Wise Error Rate p-value): The Family-Wise Error Rate (FWER) p-value is another measure of statistical significance that controls for multiple testing. It represents the probability of obtaining the observed NES or a more extreme value under the assumption of no association between the gene set and the phenotype. FWER p-values are typically more conservative than FDR q-values.

# Tag % (Tag percentage): The Tag percentage represents the percentage of genes in the input list (ranked list of all human genes) that overlap with the genes in the gene set of interest (e.g., kidney GO terms). It indicates the extent to which the gene set is represented in your input data.

# Gene % (Gene percentage): The Gene percentage represents the percentage of genes in the gene set of interest (e.g., kidney GO terms) that are found in the input list. It indicates how well the genes from your gene set are represented in your input data.


# In[13]:


# Plots
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


# In[15]:


# Plots
terms = pre_res.res2d.Term
axs = pre_res.plot(terms=terms[0])


# In[ ]:





# In[17]:


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
pre_res_kid_GO.res2d


# In[18]:


# Plot result 
terms = pre_res_kid_GO.res2d.Term

axs = pre_res_kid_GO.plot(terms=terms[0],
                   #legend_kws={'loc': (1.2, 0)}, # set the legend loc
                   show_ranking=True, 
                   figsize=(15,10)
                  )


# In[22]:


## Enrichment analysis for kidney GO terms in GO terms ranked by NES of previous analysis

# chatGPT => Problem: only checks if there are more elements of the unrankedSublist ind the upper
# n = len(unrankedSublist) than by chance. Exact ranking not relevant
# for len(kidney_GO_terms) = 

def test_statistic(ranked_list, unranked_sublist, occur_list):
    # consider top portion of ranked list
    top_portion = ranked_list[:len(unranked_sublist)]  
    
    # count occurences of elements of unrankedSublist in the top portion of the ranked list
    occurrences = sum(1 for element in unranked_sublist if element in top_portion)
    
    # append occurences list
    occur_list.append(occurrences)
    return occurrences, occur_list

def permute(ranked_list):
    # create a copy of the ranked list
    permuted_ranked_list = ranked_list[:] 
    
    # shuffle the copied list
    random.shuffle(permuted_ranked_list)  
    return permuted_ranked_list

def permutation_test(ranked_list, unranked_sublist, num_permutations, occur_list):
    random.seed(1)
    observed_statistic, occur_list = test_statistic(ranked_list, unranked_sublist, occur_list)
    permuted_statistics = []
    
    for _ in range(num_permutations):
        permuted_ranked_list = permute(ranked_list)  # Randomly permute ranked list
        permuted_statistic, occur_list = test_statistic(permuted_ranked_list, unranked_sublist, occur_list)
        permuted_statistics.append(permuted_statistic)
    
    # calculate p-value
    p_value = sum(1 for stat in permuted_statistics if stat >= observed_statistic) / num_permutations
    return p_value, occur_list



ranked_list = NGS_pre_res.sort_values(by='NES', ascending=False)['GO_term'].copy().values
unranked_sublist = [i for i in kidney_GO_terms if i in ranked_list]


p_value, occur_list = permutation_test(ranked_list = ranked_list,
                                       unranked_sublist = unranked_sublist,
                                       num_permutations = 10000, 
                                       occur_list=[]
                                      )

print("p-value:", p_value)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# NOTES: TO BE DELETED

# file_path = "/fast/users/rankn_c/work/miniconda/envs/ncs_env/lib/python3.11/site-packages/gseapy/gsea.py"

# with open(file_path, "r") as file:
#     lines = file.readlines()
#     for i, line in enumerate(lines[317:470], start=317):
#         print(f"Line {i}: {line.strip()}")



# get kidney GO terms
def get_all_children(GO_term, all_children=None, checked_terms=None):
    """Function that returns a list of all descendants of a GO term and the GO term itself"""
    
    if checked_terms is None:
        checked_terms = []
    
    if all_children is None:
        all_children = []

    all_children.append(GO_term)
        
    direct_children = get_direct_children(GO_term)
    checked_terms.append(GO_term)
    
    for child in direct_children:
        if child not in checked_terms:
        
            if get_direct_children(child) == [child]:
                all_children.append(child)
                checked_terms.append(child)

            else:
                get_all_children(child, all_children, checked_terms)        
        
    return list(set(all_children))
    
    
# get_all_children("GO:0097744")
# a = get_all_children("GO:0035812")
# a = get_all_children("GO:0097254")
a = get_all_children("GO:0072001")




def get_all_children(GO_term, all_children=None, checked_terms=None):
#     print(GO_term)
#     print(get_direct_children(GO_term))
#     print("______________")
    
    if checked_terms is None:
        checked_terms = []
    
    if all_children is None:
        all_children = []

    all_children.append(GO_term)
        
    direct_children = get_direct_children(GO_term)
    checked_terms.append(GO_term)
    
    for child in direct_children:
        if child not in checked_terms:
        
            if get_direct_children(child) == [child]:
                all_children.append(child)
                checked_terms.append(child)
                print(f"leave_node: {child}")
                print("______________")
            else:
#                 print(child)
#                 print(get_direct_children(child))
#                 print("______________")


                get_all_children(child, all_children, checked_terms)        
        
    return list(set(all_children))
    
    
# get_all_children("GO:0097744")
# a = get_all_children("GO:0035812")
# a = get_all_children("GO:0097254")
a = get_all_children("GO:0072001")

a


all_children = []

def get_all_children(GO_term, all_children=None):
    print(GO_term)
    print(get_direct_children(GO_term))
    print("______________")
    
    if all_children is None:
        all_children = []

    all_children.append(GO_term)
        
    direct_children = get_direct_children(GO_term)
    
    for child in direct_children:
        
        if get_direct_children(child) == [child]:
            all_children.append(child)
            print(f"leave_node: {child}")
            print("______________")
        else:
            print(child)
            print(get_direct_children(child))
            print("______________")
            
            
            get_all_children(child, all_children)        
        
    return list(set(all_children))
    
    
# get_all_children("GO:0097744")
# a = get_all_children("GO:0035812")
a = get_all_children("GO:0097254")
# a = get_all_children("GO:0072001")

a



# In[ ]:


# https://wiki.geneontology.org/Kidney_-_Agenda_for_Kidney_Development_Ontology_Content_Meeting#Ontology_Editing
# https://www.ebi.ac.uk/QuickGO/term/GO:0072001


# In[ ]:


gene_list = ['SYNPO2L', 'SCARA3', 'LOC100044683', 'CMBL', 'CLIC6', 'TACSTD2', 'DKKL1',
                'CSF1', 'IL13RA1', 'CITED1', 'PKD1', 'GANAB', 'PKD2', 'ALG6']

gene_sets = ['KEGG_2021_Human', 'GO_Molecular_Function_2023']


# Overrepresentation analysis
# if you are only intrested in dataframe that enrichr returned, please set outdir=None
enr = gp.enrichr(gene_list=gene_list, # or "./tests/data/gene_list.txt",
                 gene_sets=gene_sets,
                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )

enr.results.head(5)


# In[ ]:





# In[14]:


import pandas as pd
ID = 97
results_date = "2024-03-22"


NGS = pd.read_csv(f"results/NGS_predictions_ID{ID}_all_{results_date}.csv.gz")
NGS_ranked = NGS.sort_values('NGS', ascending = False).reset_index(drop = True)[['symbol', 'NGS']]
NGS_ranked = NGS_ranked.rename(columns={'symbol':'Gene', 'NGS':'Rank'}).drop_duplicates()

NGS_ranked


# In[16]:


['FOXC1', 'NPNT', 'SIX1', 'FOXJ1', 'OXSR1', 'TCF21', 'PKD2', 'BMP7', 'ROBO2', 'BMP4', 'SALL1', 'WNT11', 'WT1', 'LHX1', 'SLIT2', 'GPC3', 'CER1', 'FGFR2']


# In[97]:


# check enrichment score calculation in gseapy: focus on weight/correl_vec

import numpy as np
# gene_list = NGS_ranked.query("Rank > 0.998")['Gene'].values
gene_list = NGS_ranked['Gene'].values

N = len(gene_list) # normally 19255, here after query 15 for better illustration

# gene_set = ['PKD1', 'GATA3', 'TNFRSF11B', 'AKAP11'] # order shouldn't matter?
# gene_set = ['PKD1', 'GATA3', 'TNFRSF11B', 'PLG'] # order shouldn't matter?
gene_set = gp.get_library("GO_Biological_Process_2023")['Ureteric Bud Development (GO:0001657)']
print(gene_set)

weight = 1 # default = 1
# correl_vector = NGS_ranked.query("Rank > 0.998")['Rank'].values
correl_vector = NGS_ranked['Rank'].values 

nperm = 4 # no. permutations
seed = 2
scale = False # default
single = False # default

tag_indicator = np.in1d(gene_list, gene_set, assume_unique=True).astype(int)  # array([1, 0, 1, 1, ..., 0, 0, 0])

if weight == 0:
    correl_vector = np.repeat(1, N) # array([1, 1, 1, ..., 1, 1, 1])
else:
    correl_vector = np.abs(correl_vector) ** weight # TODO: is correl_vector rank_value??, 


hit_ind = np.flatnonzero(tag_indicator).tolist() # [0, 2, 3, 733]
hit_ind

axis = 1
tag_indicator = np.tile(tag_indicator, (nperm + 1, 1)) # shape = (nperm+1, 19255), each of the nperm elements: [1, 0, 1, 1, ..., 0, 0, 0]
correl_vector = np.tile(correl_vector, (nperm + 1, 1)) # dhape = (nperm+1, 19255), for weight = 0: each of the nperm elements: [1, 1, 1, 1, ..., 1, 1, 1]


rs = np.random.RandomState(seed)

for i in range(nperm): # only shuffles first 4 elements of tag_indicator, the last one is the "observed" one
    rs.shuffle(tag_indicator[i])
    
tag_indicator # randomly shuffled values in each of the nperm arrays in tag_indicator, not the last one

Nhint = tag_indicator.sum(axis=axis, keepdims=True) # shape(nperm+1, 1) value here 3 for all 

sum_correl_tag = np.sum(correl_vector * tag_indicator, axis=axis, keepdims=True) # shape(nperm+1, 1) value for weight=0 is here 3 for all 

no_tag_indicator = 1 - tag_indicator # shape as tag_indicator, has 0 and 1 reverse to tag_indicator

Nmiss = N - Nhint # shape = (nperm+1, 1) values here all 12 (= 15-3 or 19255-3 if NGS_ranked not queried)
norm_tag = 1.0 / sum_correl_tag # shape = (nperm+1, 1) values here all 0.33 (= 1/3)
norm_no_tag = 1.0 / Nmiss # shape = (nperm+1, 1) values here all 0.083 (= 1/12)

RES = np.cumsum(
    tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag,
    axis=axis,
) # shape = 5x15

if scale:
    RES = RES / N
if single:
    es_vec = RES.sum(axis=axis) # len = 5
else:
    max_ES, min_ES = RES.max(axis=axis), RES.min(axis=axis) # shapes = 5x1
    es_vec = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES) # shape = 5x1
    # extract values
es, esnull, RES = es_vec[-1], es_vec[:-1], RES[-1, :]




# In[ ]:


def enrichment_score(
        self,
        gene_list: Iterable[str],
        correl_vector: Iterable[float],
        gene_set: Dict[str, List[str]],
        weight: float = 1.0,
        nperm: int = 1000,
        seed: int = 123,
        single: bool = False,
        scale: bool = False,
    ):
        """This is the most important function of GSEApy. It has the same algorithm with GSEA and ssGSEA.

        :param gene_list:       The ordered gene list gene_name_list, rank_metric.index.values
        :param gene_set:        gene_sets in gmt file, please use gmt_parser to get gene_set.
        :param weight:  It's the same with gsea's weighted_score method. Weighting by the correlation
                                is a very reasonable choice that allows significant gene sets with less than perfect coherence.
                                options: 0(classic),1,1.5,2. default:1. if one is interested in penalizing sets for lack of
                                coherence or to discover sets with any type of nonrandom distribution of tags, a value p < 1
                                might be appropriate. On the other hand, if one uses sets with large number of genes and only
                                a small subset of those is expected to be coherent, then one could consider using p > 1.
                                Our recommendation is to use p = 1 and use other settings only if you are very experienced
                                with the method and its behavior.

        :param correl_vector:   A vector with the correlations (e.g. signal to noise scores) corresponding to the genes in
                                the gene list. Or rankings, rank_metric.values
        :param nperm:           Only use this parameter when computing esnull for statistical testing. Set the esnull value
                                equal to the permutation number.
        :param seed:            Random state for initializing gene list shuffling. Default: seed=None

        :return:

        ES: Enrichment score (real number between -1 and +1)

        ESNULL: Enrichment score calculated from random permutations.

        Hits_Indices: Index of a gene in gene_list, if gene is included in gene_set.

        RES: Numerical vector containing the running enrichment score for all locations in the gene list .

        """
        N = len(gene_list)
        # Test whether each element of a 1-D array is also present in a second array
        # It's more intuitive here than original enrichment_score source code.
        # use .astype to covert bool to integer
        tag_indicator = np.in1d(gene_list, gene_set, assume_unique=True).astype(
            int
        )  # notice that the sign is 0 (no tag) or 1 (tag)

        if weight == 0:
            correl_vector = np.repeat(1, N)
        else:
            correl_vector = np.abs(correl_vector) ** weight

        # get indices of tag_indicator
        hit_ind = np.flatnonzero(tag_indicator).tolist()
        # if used for compute esnull, set esnull equal to permutation number, e.g. 1000
        # else just compute enrichment scores
        # set axis to 1, because we have 2D array
        axis = 1
        tag_indicator = np.tile(tag_indicator, (nperm + 1, 1))
        correl_vector = np.tile(correl_vector, (nperm + 1, 1))
        # gene list permutation
        rs = np.random.RandomState(seed)
        for i in range(nperm):
            rs.shuffle(tag_indicator[i])
        # np.apply_along_axis(rs.shuffle, 1, tag_indicator)

        Nhint = tag_indicator.sum(axis=axis, keepdims=True)
        sum_correl_tag = np.sum(correl_vector * tag_indicator, axis=axis, keepdims=True)
        # compute ES score, the code below is identical to gsea enrichment_score method.
        no_tag_indicator = 1 - tag_indicator
        Nmiss = N - Nhint
        norm_tag = 1.0 / sum_correl_tag
        norm_no_tag = 1.0 / Nmiss

        RES = np.cumsum(
            tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag,
            axis=axis,
        )

        if scale:
            RES = RES / N
        if single:
            es_vec = RES.sum(axis=axis)
        else:
            max_ES, min_ES = RES.max(axis=axis), RES.min(axis=axis)
            es_vec = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES)
        # extract values
        es, esnull, RES = es_vec[-1], es_vec[:-1], RES[-1, :]

        return es, esnull, hit_ind, RES

