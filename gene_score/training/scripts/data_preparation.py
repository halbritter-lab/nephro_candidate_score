#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.environ['CONFIG_FILE'] = '/fast/work/users/rankn_c/halbritter/nephro_candidate_score/config_NCS.yml' 
# TODO: change


# In[2]:


# import basic modules
import numpy as np
import pandas as pd
import os
import sys
import yaml


# import preprocessing modules
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA


# In[3]:


# get config file
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

# append path where common functions are located
sys.path.append(f"{config_vars['ML_projectsdir']}{project_name}")

# import common functions
from common_functions.training_helper_functions import *

# define score type
score = 'gs'
score_string, id_string, label_string = get_score_specific_args(score)


# In[ ]:





# In[11]:


# load features and labels
raw_feat = pd.read_csv(f"gene_score/features/results/gene_features_{config_vars['creation_date_gs']}.csv.gz")
pos_genes = pd.read_csv(f"gene_score/labels/results/positive_genes_{config_vars['creation_date_gs']}.csv.gz")
neg_genes = pd.read_csv(f"gene_score/labels/results/dispensible_genes_{config_vars['creation_date_gs']}.csv.gz")

# load table with filling methods for missing values
filling_methods = pd.read_csv(f"gene_score/training/raw/fill_missing_methods_{config_vars['creation_date_gs']}.csv.gz", delimiter=";")


# In[17]:


# create new column indicating if a evidence_count > 1
pos_genes['ec2345'] = np.where(pos_genes['evidence_count'] > 1, 1, 0)

# filter only positive genes with evidence_count > 1
pos_genes = pos_genes.query("ec2345 == 1")

# set negative genes to ec2345 = 0 
neg_genes['ec2345'] = 0

# combine positive and negative genes
labels = pd.concat([neg_genes[['hgnc_id_int', 'ec2345']], pos_genes[['hgnc_id_int', 'ec2345']]])


# In[24]:


# split labels into 80% training and 20% test set stratified by ec2345
train_ids, test_ids, _, _ = train_test_split(labels.drop(columns=['ec2345']), labels['ec2345'], test_size=0.2, stratify=labels['ec2345'], random_state=42)

train_ids_arr = train_ids['hgnc_id_int'].values
test_ids_arr = test_ids['hgnc_id_int'].values

# get training labels df
labels_train = labels[labels['hgnc_id_int'].isin(train_ids_arr)]

# get training features df
feat_train = raw_feat[raw_feat['hgnc_id_int'].isin(train_ids_arr)]

# get test labels df
labels_test = labels[labels['hgnc_id_int'].isin(test_ids_arr)]

# get test features df
feat_test = raw_feat[raw_feat['hgnc_id_int'].isin(test_ids_arr)]

# save training and test data
labels_train.to_csv(f"gene_score/training/train_test_data/labels_train_{config_vars['data_prep_date_gs']}.csv.gz", index=False, compression='gzip')
labels_test.to_csv(f"gene_score/training/train_test_data/labels_test_{config_vars['data_prep_date_gs']}.csv.gz", index=False, compression='gzip')
feat_train.to_csv(f"gene_score/training/train_test_data/feat_train_{config_vars['data_prep_date_gs']}.csv.gz", index=False, compression='gzip')
feat_test.to_csv(f"gene_score/training/train_test_data/feat_test_{config_vars['data_prep_date_gs']}.csv.gz", index=False, compression='gzip')

# create table with median values of training set for filling NA
median_values_train = feat_train.drop(columns = ['hgnc_id_int']).median().reset_index()
median_values_train.columns = ['feature', 'median']

# save median df
median_values_train.to_csv(f"gene_score/training/train_test_data/median_values_train_{config_vars['data_prep_date_gs']}.csv.gz", index=False, compression='gzip')


# In[ ]:


# Note: labels_train contains more genes than feat_train as also some pseudogenes are included in labels_train. 
# These are filtered out in raw_feat (and so in feat_train). In the training process later, labels_train and
# feat_train are merged by an inner join joined on 'hgnc_id_int' (inner join), so no problem arises.


# In[ ]:





# In[14]:


## Create a training set that contains the most important principal components of each feature group, that 
# has highly correlated features

# all feature groups
feature_groups = ['gnomad', 'cellxgene', 'descartes', 'gtex', 'mgi', 'paralogues', 'phasCons', 'CpG_o2e']

# feature groups on which PCA should be performed
groups_for_PCA = [i for i in feature_groups if i not in ['mgi', 'paralogues', 'phasCons', 'CpG_o2e']] 
   
# create a dictionary that stores how many PC should be kept for each group to reach cumulative variance explained > 95% (99%)    
n_comp_dic_95 = {}
n_comp_dic_99 = {}

# manually analyse PCA for each group and select number of components to keep
for group in groups_for_PCA:
    reduced_df, loadings_df, cumsum = PCA_on_feature_group(
        group=group,
        model='logReg',
        plot=False,
        n_components=False,
        omit_scaling_groups=['mgi', 'paralogues'],
        score=score
    )  

    # keep number of components that reach cumulative variance explained > 95%
    n_comp_dic_95[group] = sum([i < 0.95 for i in cumsum]) + 1

    # keep number of components that reach cumulative variance explained > 99%
    n_comp_dic_99[group] = sum([i < 0.99 for i in cumsum]) + 1


# build a new training feature df with relevant principal components (95% variance kept/99% variance kept)
# and all features on which PCA was not performed (e.g. MGI, paralogues)
feat_train = pd.read_csv(f"{score_string}/training/train_test_data/feat_train_{config_vars[f'data_prep_date_{score}']}.csv.gz")
feat_train_reduced_95 = feat_train[['hgnc_id_int'] + get_features_from_groups(['mgi', 'paralogues', 'phasCons', 'CpG_o2e'], feat_train, score)]
feat_train_reduced_99 = feat_train[['hgnc_id_int'] + get_features_from_groups(['mgi', 'paralogues', 'phasCons', 'CpG_o2e'], feat_train, score)]


for group in groups_for_PCA:
    reduced_df_95, loadings_df_95, cumsum = PCA_on_feature_group(
        group=group,
        model='logReg',
        plot=False,
        n_components=n_comp_dic_95[group],
        omit_scaling_groups=['mgi', 'paralogues'],
        score=score
    )  
    
    # merge with original features from groups that were not reduced with PCA
    feat_train_reduced_95 = pd.merge(feat_train_reduced_95, reduced_df_95, on='hgnc_id_int')
    feat_train_reduced_95.to_csv(f"gene_score/training/train_test_data/feat_train_reduced_95_{data_prep_date_gs}.csv", index=False)

    
    reduced_df_99, loadings_df, cumsum = PCA_on_feature_group(
        group=group,
        model='logReg',
        plot=False,
        n_components=n_comp_dic_99[group],
        omit_scaling_groups=['mgi', 'paralogues'],
        score=score
    )  
    
    # merge with original features from groups that were not reduced with PCA
    feat_train_reduced_99 = pd.merge(feat_train_reduced_99, reduced_df_99, on='hgnc_id_int')
    feat_train_reduced_99.to_csv(f"gene_score/training/train_test_data/feat_train_reduced_99_{data_prep_date_gs}.csv", index=False)
  
    
    


# In[ ]:





# In[ ]:




