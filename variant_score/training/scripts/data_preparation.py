#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.environ['CONFIG_FILE'] = '/fast/work/users/rankn_c/halbritter/nephro_candidate_score/gene_score/training/config_NCS.yml' # TODO: change
# TODO: change

# in terminal:
# export CONFIG_FILE=/fast/work/users/rankn_c/halbritter/nephro_candidate_score/gene_score/training/config_NCS.yml


# In[2]:


# import basic modules
import numpy as np
import pandas as pd
import yaml
import os
import vcf
import gzip
import re
from datetime import datetime

# import preprocessing modules
from sklearn.model_selection import train_test_split

# set options
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

# load config file
CONFIG_FILE = os.getenv('CONFIG_FILE')

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"
script_path = "/variant_score/"

# read configs
with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}{script_path}")


# In[5]:


# load features and labels
features_labels = pd.read_csv(f"features_labels/results/features_labels_{config_vars['data_prep_date_vs']}.csv.gz")

# split labels into 80% training and 20% test set stratified by P_LP
labels_train, labels_test, _, _ = train_test_split(features_labels[['var_ID', 'P_LP']], features_labels['P_LP'], test_size=0.2, stratify=features_labels['P_LP'], random_state=42)

labels_train = labels_train.reset_index(drop=True)
labels_test = labels_test.reset_index(drop=True)


# get training features
feat_train = pd.merge(labels_train[['var_ID']], features_labels, how='left', left_on='var_ID', right_on='var_ID').drop(columns=['P_LP'])

# get testing features
feat_test = pd.merge(labels_test[['var_ID']], features_labels, how='left', left_on='var_ID', right_on='var_ID').drop(columns=['P_LP'])


# save training and test data
labels_train.to_csv(f"training/train_test_data/labels_train_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
labels_test.to_csv(f"training/train_test_data/labels_test_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
feat_train.to_csv(f"training/train_test_data/feat_train_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
feat_test.to_csv(f"training/train_test_data/feat_test_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')

# create table with median values of training set for filling NA
median_values_train = feat_train.drop(columns=['var_ID']).median().reset_index()
median_values_train.columns = ['feature', 'median']

# save median df
median_values_train.to_csv(f"training/train_test_data/median_values_train_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')


