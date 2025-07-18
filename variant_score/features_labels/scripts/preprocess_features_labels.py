#!/usr/bin/env python
# coding: utf-8

#### Preprocess features and labels ####
# This script loads the raw VEP-annotated ClinVar variants, analyses the features (e.g. check for missing values),
# preprocess the features and outputs the preprocessed final features_labels file.
# converted from .ipynb

# import basic modules
import numpy as np
import pandas as pd
import yaml
import urllib.request
import os
import vcf
import gzip
import re
from datetime import datetime
import os


# set options
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)


# CONFIG_FILE = "config_NCS.yml"
CONFIG_FILE = os.getenv('CONFIG_FILE')

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"

# read configs
with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}")


# In[3]:


# load extracted features from VEP annotation

# create empty df
raw_feat = pd.DataFrame()

# concatenate features for all chromosomes (exclude MT variants for now)
all_chrom = list(np.arange(22) + 1) + ['X', 'Y']


for chrom in all_chrom:
    print(chrom)
    chrom_feat = pd.read_csv(f"variant_score/features_labels/results/raw_features_clinvar_vars_kid-gen_2345_chr{chrom}_{config_vars['raw_features_clinvar_kid-gen_date_vs']}.csv.gz", low_memory=False)
    
    # calculate maximum length of REF and ALT for each row and remove REF or ALT values longer than 50 characters
    if chrom_feat.shape[0] > 0:
        chrom_feat['max_length_REF_ALT'] = chrom_feat.apply(lambda row: max(len(row['REF']), len(row['ALT'])), axis=1)
        chrom_feat = chrom_feat.query("max_length_REF_ALT <= 50").drop(columns=['max_length_REF_ALT'])
    
    raw_feat = pd.concat([raw_feat, chrom_feat], axis=0).reset_index(drop=True)


# create an own variant ID for each record: CHROM_POS_REF_ALT
raw_feat = raw_feat.copy()
raw_feat.loc[:, 'var_ID'] = raw_feat['CHROM'].astype(str) + '_' \
+ raw_feat['POS'].astype(str) + '_' + raw_feat['REF'].astype(str) \
+ '_' + raw_feat['ALT'].astype(str)


# map IMPACT to integer
impact_mapping = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
raw_feat['IMPACT_num'] = raw_feat['IMPACT'].map(impact_mapping)


# In[4]:


# # check for which features there > 1 value per variant
# unequal_entries = []

# len_unique_var_IDs = len(raw_feat['var_ID'].unique())
# for col in raw_feat.columns:

#     column_to_check = col
#     # raw_feat[['var_ID', column_to_check]].drop_duplicates()
#     a = raw_feat[['var_ID', column_to_check]].drop_duplicates().shape[0]

#     print(column_to_check)
#     print(a)
#     print(len_unique_var_IDs)
#     if a != len_unique_var_IDs:
#         unequal_entries = unequal_entries + [col]
        
#     print("-------------")
    
# unequal_entries


# In[5]:


# # manually check variants with multiple values per same variant for different features
# column_to_check = 'phastCons100way_vertebrate_rankscore'
# df = raw_feat[['var_ID', column_to_check]].drop_duplicates()
# df.groupby('var_ID').size().reset_index(name='count').query("count > 1")


# In[6]:


# count_gr1 = df.query(f'{column_to_check}.notna()').groupby('var_ID').size().reset_index(name='count').query("count > 1")
# count_gr1


# In[7]:


# # get pairs of genes with same variant
# gene_pairs = []

# count = 0

# for i in count_gr1['var_ID'].unique():
#     var = [i]
#     symbol_pair = list(raw_feat.query("var_ID in @var")['SYMBOL'].unique())
#     gene_pairs = gene_pairs + [symbol_pair]
    
#     count = count + 1
#     if count % 100 == 0:
#         print(list({tuple(sublist) for sublist in gene_pairs}))


# In[ ]:


# double_varIds = df.query(f'{column_to_check}.notna()').groupby('var_ID').size().reset_index(name='count').query("count > 1")['var_ID']

# raw_feat.query("var_ID == '17_19357902_C_T'")
# # raw_feat.query("var_ID in @ double_varIds")[['var_ID', 'HGNC_ID', 'SYMBOL', 'hgnc_id_int']]


# In[9]:


# get unique CLNSIG values
# raw_feat[['var_ID', 'CLNSIG']].drop_duplicates().groupby('CLNSIG').size().reset_index(name='count')


# In[ ]:





# In[ ]:





# In[10]:


## FINDINGS:
# - Consequence, IMPACT => different values for same variant
# - same variant can be annotated with 2 different SYMBOL/HGNC_ID/hgnc_id_int if two genes (both kidney-genetic 
# genes overlap), e.g. 10_70888586_G_T in PCBD1 and SGPL1 
# overlapping genes: [('EVC', 'EVC2'), ('NADSYN1', 'DHCR7'), ('ABCC8', 'KCNJ11'), ('DAAM2', 'MOCS1'), ('SGPL1', 'PCBD1'), ('PCBD1', 'SGPL1'), ('UQCC2', 'ITPR3'), ('COL4A3', 'COL4A4'), ('ITPR3', 'UQCC2'), ('TSC2', 'PKD1'), ('KCNJ11', 'ABCC8')]
# - NMD: NA or e.g. NMD_escaping_variant  for same variant
# - SpliceAI_pred_DS_...: NA or other value for same variant
# - BLOSOM62: different values for same variant, e.g. NA, 0, -2
# - MaxEntScan_alt, MaxEntScan_diff, MaxEntScan_ref: different values for same variant
# - REVEL_rankscore: different values for same variant
# - PrimateAI_rankscore: different values only in one case 17_19357902_C_T (Uncertain_significance)
# - phastCons100way_vertebrate_rankscore: NA or other value for same variant


# REVEL_rankscore: the higher the more likely pathogenic
# MaxEntScan_ref: This refers to the MaxEntScan score calculated for the reference (normal) sequence surrounding the splice site.
# MaxEntScan_alt: This refers to the MaxEntScan score calculated for the alternative (mutated) sequence surrounding the splice site. It helps assess how the mutation may affect splice site strength compared to the reference sequence.
# MaxEntScan_diff: This represents the difference between the MaxEntScan scores of the alternative and reference sequences. A positive difference indicates that the mutation strengthens the splice site, while a negative difference suggests a weakening of the splice site.


# In[11]:


# filter only for B/LB, P/LP
selected_CLNSIG = ['Benign/Likely_benign', 'Likely_benign', 'Benign', 'Pathogenic/Likely_pathogenic',
       'Likely_pathogenic', 'Pathogenic']

raw_feat = raw_feat.query("CLNSIG in @selected_CLNSIG")


# In[12]:


# check missing value proportion
col_to_check = 'MaxEntScan_ref'

not_missing = len(raw_feat.query(f'{col_to_check}.notna()')['var_ID'].unique())
not_missing_prop = not_missing/len(raw_feat['var_ID'].unique())
not_missing_prop

def get_missing_val_prop(col_to_check):
    not_missing = len(raw_feat.query(f'{col_to_check}.notna()')['var_ID'].unique())
    not_missing_prop = not_missing/len(raw_feat['var_ID'].unique())
    return(1-not_missing_prop)

for col in raw_feat.columns:
    print(f"missing value proportion for {col}: {get_missing_val_prop(col)}")


# In[ ]:





# In[14]:


# get single feature 

def get_single_feature_df(feature, method, ms_aggregation_symbol=None, ms_prefix=None):
    """
    Selects for given feature one value for each variant defined by 'method'.
    - 'method' has to be 'notna', 'max', 'min'.
    - 'ms_aggregation_symbol' is only needed if 'method' == 'multiple_split' and has to be the symbol 
    that separates multiple values per row of the given feature, e.g. "&", "|".
    - 'ms_prefix' is only needed if 'method' == 'multiple_split' and is the prefix that is added before each column of distinct values.
    Returns: pd.dataframe with columns 'varID' and {feature}
    """
    
    feat_df = raw_feat[['var_ID', feature]].drop_duplicates()
    
    if method == 'notna': # get not NaN value per variant
        feat_df = feat_df.query(f'{feature}.notna()').drop_duplicates()
        # check if there are multiple not NaN values for one variant
        if feat_df.groupby('var_ID').size().reset_index(name='count').query("count >1").shape[0] > 1:
            raise ValueError("Multiple not NaN values for some variants. Choose other method.")
        else: return(feat_df)
     
    elif method == 'max': # get maximum value per variant (NaN will be ignored if other values are available)
        return(feat_df.groupby('var_ID').max().reset_index())
    
    elif method == 'min': # get minimum value per variant (NaN will be ignored if other values are available)
        return(feat_df.groupby('var_ID').min().reset_index())
    
    elif method == 'multiple_split': # splits multiple values per variant in separate columns
        
        # if only one not NaN value exists for feature
        if len(feat_df.query(f"{feature}.notna()")[feature].unique()) == 1:
            unique_feature_val = feat_df.query(f"{feature}.notna()")[feature].unique()[0]
            
            feat_df = feat_df.query(f'{feature}.notna()').drop_duplicates()

            feat_mapping = {unique_feature_val: 1}
            feat_df[unique_feature_val] = feat_df[feature].map(feat_mapping)
            feat_df = feat_df.drop(columns=[feature])
            return(feat_df)

        # if multiple not NaN values exist for feature
        else:
            if ms_aggregation_symbol is None:
                raise ValueError("Provide 'ms_aggregation_symbol' argument.")
            # aggregate all values per variant of the given feature
            feat_df = feat_df.groupby('var_ID')[feature].agg(f"{ms_aggregation_symbol}".join).reset_index()

            # split the feature column into multiple separate columns for each distinct value, respectively
            dummy_columns = feat_df[feature].str.get_dummies(sep=ms_aggregation_symbol)
            if ms_prefix is None:
                ms_prefix = feature
            dummy_columns = dummy_columns.add_prefix(f"{ms_prefix}_")

            # concatenate the dummy columns with the original DataFrame
            return(pd.concat([feat_df, dummy_columns], axis=1))
    
    else:
        raise ValueError("'method' argument must be 'notna', 'max', 'min', or 'multiple_split'.")
        
# a = get_single_feature_df(feature='Consequence', 
#                       method='multiple_split', 
#                       ms_aggregation_symbol = "&",
#                       ms_prefix="Csq"
#                      )   



# In[15]:


# import pandas as pd

# feature = 'NMD'

# # Assuming your DataFrame is named df
# feat_df = raw_feat[['var_ID', feature]].drop_duplicates()


# unique_feature_val = feat_df.query(f"{feature}.notna()")[feature].unique()[0]

# # Create dummy variables for the 'NMD' column
# # dummy_df = pd.get_dummies(feat_df['NMD'])

# # Select the 'NMD_escaping_variant' column
# # unique_val_column = dummy_df[unique_feature_val]

# feat_mapping = {unique_feature_val: 1}
# feat_df[unique_feature_val] = feat_df[feature].map(feat_mapping)



# # Rename the column to 'NMD_escaping_variant' and assign it back to the DataFrame
# # feat_df[unique_feature_val] = unique_val_column

# # Fill NaN values with NaN (optional, but may be useful if there are other columns with NaN values)
# # feat_df[feat_df].fillna(value=pd.NA, inplace=True)

# # Drop the original 'NMD' column if necessary
# # df.drop(columns=['NMD'], inplace=True)

# # feat_df.query("NMD.notna()")

# # df['my_column'] = df['my_column'].astype(int)
# feat_df.query("NMD.notna()")


# In[17]:


features_labels = raw_feat[['var_ID', 'CLNSIG']].drop_duplicates()

## map CLNSIG values to 0 and 1

def map_pathogenicity(sig):
    if sig in ['Benign/Likely_benign', 'Likely_benign', 'Benign']:
        return 0
    elif sig in ['Pathogenic/Likely_pathogenic', 'Likely_pathogenic', 'Pathogenic']:
        return 1
    else:
        raise ValueError(f"No mapping for value '{sig}'. Adapt map_pathogenicity().")

# apply the function to create the 'P_LP' column
features_labels['P_LP'] = raw_feat['CLNSIG'].apply(lambda x: map_pathogenicity(x))


# define features that should be included for training of the variant score
selected_features = ['Consequence', 
                     'gnomADe_AF',
                     'gnomADg_AF', 
                     'CADD_PHRED', 
                     'NMD',
                     'SpliceAI_pred_DS_AG', 
                     'SpliceAI_pred_DS_AL', 
                     'SpliceAI_pred_DS_DG',
                     'SpliceAI_pred_DS_DL',
                     'phastCons100way_vertebrate_rankscore',
                     'IMPACT_num']

# define methods/arguments to choose a single value for each variant if multiple values are available for a variant
# [method, ms_aggregation_symbol, ms_prefix]
feature_args = {
    'Consequence': ['multiple_split', '&', 'Csq'],
    'IMPACT_num': ['max', None, None],
    'gnomADe_AF': ['notna', None, None],
    'gnomADg_AF': ['notna', None, None],
    'CADD_PHRED': ['notna',  None, None],
    'NMD': ['multiple_split',  None, None],
    'SpliceAI_pred_DS_AG': ['notna', None, None],
    'SpliceAI_pred_DS_AL': ['notna', None, None],
    'SpliceAI_pred_DS_DG': ['notna', None, None],
    'SpliceAI_pred_DS_DL': ['notna', None, None],
    'phastCons100way_vertebrate_rankscore': ['notna', None, None]
}


# add each feature with the defined method to chose from multiple values per variant
for feat in selected_features:
    single_feat = get_single_feature_df(feature=feat, 
                      method=feature_args[feat][0], 
                      ms_aggregation_symbol = feature_args[feat][1],
                      ms_prefix=feature_args[feat][2]
                     ) 
    
    
    features_labels = pd.merge(features_labels, single_feat, how='left', left_on='var_ID', right_on='var_ID')
    


# In[20]:


# remove variants that have only missing values in features
features = list(set(features_labels.columns) - set(['var_ID', 'CLNSIG', 'P_LP']))
features_labels = features_labels.dropna(subset=features, how='all')

# drop unnecessary columns
features_labels = features_labels.drop(columns=['CLNSIG', 'Consequence'])

# write csv
current_date = datetime.today().strftime('%Y-%m-%d')
# current_date = "2024-04-08"
features_labels.to_csv(f"variant_score/features_labels/results/features_labels_{current_date}.csv.gz", index=False, compression='gzip')



# In[ ]:




