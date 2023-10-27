#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import config
from config_ML import *

# import basic modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import preprocessing modules
from sklearn.model_selection import train_test_split


# set options
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', 500)


# In[2]:


# load features and labels
raw_feat = pd.read_csv(f'{raw_data_dir}/gene_features_{creation_date}.csv')
pos_genes = pd.read_csv(f'{raw_data_dir}/positive_genes_{creation_date}.csv')
neg_genes = pd.read_csv(f'{raw_data_dir}/dispensible_genes_{creation_date}.csv')

# load table with filling methods for missing values
filling_methods = pd.read_csv(f'{raw_data_dir}/fill_missing_methods_{creation_date}.csv', delimiter=";")


# In[3]:


#### LABELS ####

# create new column indicating if a evidence_count > 1
pos_genes['ec2345'] = np.where(pos_genes['evidence_count'] > 1, 1, 0)

# filter only positive genes with evidence_count > 1
pos_genes = pos_genes.query("ec2345 == 1")

# set negative genes to ec2345 = 0 
neg_genes['ec2345'] = 0

# combine positive and negative genes
labels = pd.concat([neg_genes[['hgnc_id_int', 'ec2345']], pos_genes[['hgnc_id_int', 'ec2345']]])

# print out number of instances
print("Instances with ec2345 == 1: ", len(labels[labels['ec2345'] == 1]))
print("Instances with ec2345 == 0: ", len(labels[labels['ec2345'] == 0]))
print(labels.shape)


# In[4]:


# split labels into 80% training and 20% test set stratified by ec2345
train_ids, test_ids, _, _ = train_test_split(labels.drop(columns=['ec2345']), labels['ec2345'], test_size=0.2, stratify=labels['ec2345'], random_state=42)

print('shape train_ids: ', train_ids.shape)
print('shape test_ids: ', test_ids.shape)

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

print('shape labels_test: ', labels_test.shape)
print('shape labels_train: ', labels_train.shape)
print('shape feat_test: ', feat_test.shape)
print('shape feat_train: ', feat_train.shape)

# save training and test data
labels_train.to_csv(f'{train_test_data_dir}/labels_train_{data_prep_date}.csv', index=False)
labels_test.to_csv(f'{train_test_data_dir}/labels_test_{data_prep_date}.csv', index=False)
feat_train.to_csv(f'{train_test_data_dir}/feat_train_{data_prep_date}.csv', index=False)
feat_test.to_csv(f'{train_test_data_dir}/feat_test_{data_prep_date}.csv', index=False)


## create table with median values of training set for filling NA
# calculate the median for each column and create a new df
median_values_train = feat_train.drop(columns = ['hgnc_id_int']).median().reset_index()
median_values_train.columns = ['feature', 'median']

# save median df
median_values_train.to_csv(f'{train_test_data_dir}/median_values_train_{data_prep_date}.csv', index=False)

