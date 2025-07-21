#!/usr/bin/env python
# coding: utf-8

#### Data preparation for training of N-VS model ####
# creates 80%/20% train/test set, creates train sets with balanced IMPACT, original proportion
# converted from .ipynb

# import basic modules
import numpy as np
import pandas as pd
import yaml
import os
import vcf
import gzip
import re
from datetime import datetime
import matplotlib.pyplot as plt

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

# read configs
with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}")


# In[3]:


# load features and labels
features_labels = pd.read_csv(f"variant_score/features_labels/results/features_labels_{config_vars['data_prep_date_vs']}.csv.gz")

# split labels into 80% training and 20% test set stratified by P_LP
labels_train, labels_test, _, _ = train_test_split(features_labels[['var_ID', 'P_LP']], features_labels['P_LP'], test_size=0.2, stratify=features_labels['P_LP'], random_state=42)

labels_train = labels_train.reset_index(drop=True)
labels_test = labels_test.reset_index(drop=True)


# get training features
feat_train = pd.merge(labels_train[['var_ID']], features_labels, how='left', left_on='var_ID', right_on='var_ID').drop(columns=['P_LP'])

# get testing features
feat_test = pd.merge(labels_test[['var_ID']], features_labels, how='left', left_on='var_ID', right_on='var_ID').drop(columns=['P_LP'])


# # save training and test data
# labels_train.to_csv(f"variant_score/training/train_test_data/labels_train_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
# labels_test.to_csv(f"variant_score/training/train_test_data/labels_test_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
# feat_train.to_csv(f"variant_score/training/train_test_data/feat_train_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
# feat_test.to_csv(f"variant_score/training/train_test_data/feat_test_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')

# create table with median values of training set for filling NA
median_values_train = feat_train.drop(columns=['var_ID']).median().reset_index()
median_values_train.columns = ['feature', 'median']

# # save median df
# median_values_train.to_csv(f"variant_score/training/train_test_data/median_values_train_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')



# In[4]:


plt.hist(feat_train.IMPACT_num) 


# In[5]:


plt.hist(labels_train.P_LP) 


# In[6]:


# build a training and test set that with the same number of instances per IMPACT value
# load features and labels
features_labels = pd.read_csv(f"variant_score/features_labels/results/features_labels_{config_vars['data_prep_date_vs']}.csv.gz")

# group the DataFrame by 'IMPACT_num'
features_labels_grouped = features_labels.groupby('IMPACT_num')

# find the minimum number of rows for any group
min_rows = features_labels_grouped.size().min()

# sample the same number of rows for each IMPACT group
features_labels_sub = features_labels_grouped.apply(lambda x: x.sample(n=min_rows)).reset_index(drop=True)
 
# create a combined columen for IMPACT and P_LP for stratification
features_labels_sub['IMPACT_P_LP'] = features_labels_sub['IMPACT_num'].astype(str) + "_" + features_labels_sub['P_LP'].astype(str)

# split labels into 80% training and 20% test set stratified by P_LP
labels_train_sub, labels_test_sub, _, _ = train_test_split(features_labels_sub[['var_ID', 'P_LP']], features_labels_sub['P_LP'], test_size=0.2, stratify=features_labels_sub['IMPACT_P_LP'], random_state=42)


labels_train_sub = labels_train_sub.reset_index(drop=True)
labels_test_sub = labels_test_sub.reset_index(drop=True)


# get training features
feat_train_sub = pd.merge(labels_train_sub[['var_ID']], features_labels_sub, how='left', left_on='var_ID', right_on='var_ID').drop(columns=['P_LP', 'IMPACT_P_LP'])

# get testing features
feat_test_sub = pd.merge(labels_test_sub[['var_ID']], features_labels_sub, how='left', left_on='var_ID', right_on='var_ID').drop(columns=['P_LP', 'IMPACT_P_LP'])


# # save training and test data
# labels_train_sub.to_csv(f"variant_score/training/train_test_data/labels_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
# labels_test_sub.to_csv(f"variant_score/training/train_test_data/labels_test_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
# feat_train_sub.to_csv(f"variant_score/training/train_test_data/feat_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
# feat_test_sub.to_csv(f"variant_score/training/train_test_data/feat_test_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')

# create table with median values of training set for filling NA
median_values_train_IMPACT_prop = feat_train_sub.drop(columns=['var_ID']).median().reset_index()
median_values_train_IMPACT_prop.columns = ['feature', 'median']

# # save median df
# median_values_train_IMPACT_prop.to_csv(f"variant_score/training/train_test_data/median_values_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')


# In[7]:


print(feat_train_sub.shape)
print(labels_train_sub.shape)
print(feat_test_sub.shape)
print(labels_test_sub.shape)
print(feat_train_sub.shape[0] + feat_test_sub.shape[0])


# In[ ]:





# In[9]:


# for testing IMPACT groups separately take all remaining variants of each group

features_labels = pd.read_csv(f"variant_score/features_labels/results/features_labels_{config_vars['data_prep_date_vs']}.csv.gz")
print("features_labels: ", len(features_labels))

# load training labels 
labels_train_sub = pd.read_csv(f"variant_score/training/train_test_data/labels_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz")
labels_test_sub = pd.read_csv(f"variant_score/training/train_test_data/labels_test_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz")
feat_train_sub = pd.read_csv(f"variant_score/training/train_test_data/feat_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz")
feat_test_sub = pd.read_csv(f"variant_score/training/train_test_data/feat_test_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz")

print("labels_train_sub: ", len(labels_train_sub))
print("feat_train_sub: ", len(feat_train_sub))
print("labels_test_sub: ", len(labels_test_sub))
print("feat_test_sub: ", len(feat_test_sub))

train_var_IDs = feat_train_sub['var_ID'].values
print("train_var_IDs: ", len(train_var_IDs))

remain_features_labels = features_labels.query("var_ID not in @train_var_IDs")


for i in [1,2,3,4]:
    IMPACT_test_data = remain_features_labels.query(f"IMPACT_num == {i}")
    
    # save data
    feat_test_remain = IMPACT_test_data.drop(columns=['P_LP'])
    labels_test_remain = IMPACT_test_data[['var_ID', 'P_LP']]
    
    feat_test_remain.to_csv(f"variant_score/training/train_test_data/feat_test_remain_IMPACT_{i}_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
    labels_test_remain.to_csv(f"variant_score/training/train_test_data/labels_test_remain_IMPACT_{i}_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
    
    print(len(feat_test_remain))
    print(len(labels_test_remain))



# In[10]:


## for global AUC testing, get a test set with the original proportion of IMPACT and P_LP

# load all features and labels
features_labels = pd.read_csv(f"variant_score/features_labels/results/features_labels_{config_vars['data_prep_date_vs']}.csv.gz")
print("features_labels: ", len(features_labels))

# create a combined columen for IMPACT and P_LP for stratification
features_labels['IMPACT_P_LP'] = features_labels['IMPACT_num'].astype(str) + "_" + features_labels['P_LP'].astype(str)

# get proportion of subgroups
features_labels_grouped = features_labels.groupby('IMPACT_P_LP').size().reset_index(name='counts')

# Calculate the total number of rows
total_rows = len(features_labels)

# Calculate the proportions
features_labels_grouped['orig_prop'] = features_labels_grouped['counts'] / total_rows

features_labels_grouped


# In[11]:


## sample a test set with the same proportion of IMPACT and P_LP as in the original distribution

# in the remaining set (not used for training) create a combined column for IMPACT and P_LP for stratification
remain_features_labels['IMPACT_P_LP'] = remain_features_labels['IMPACT_num'].astype(str) + "_" + remain_features_labels['P_LP'].astype(str)

remain_features_labels_grouped = remain_features_labels.groupby('IMPACT_P_LP').size().reset_index(name='counts')

total_rows_remain = len(remain_features_labels)

remain_features_labels_grouped['prop'] = remain_features_labels_grouped['counts'] / total_rows_remain

# join with original proportion and count proportion ratio
remain_features_labels_grouped = pd.merge(remain_features_labels_grouped, features_labels_grouped[['IMPACT_P_LP', 'orig_prop']], how='inner', left_on='IMPACT_P_LP', right_on='IMPACT_P_LP')
remain_features_labels_grouped['prop_ratio'] = remain_features_labels_grouped['prop'] / remain_features_labels_grouped['orig_prop'] 

# get mininum of proportion ratio
min_prop_ratio_remain = remain_features_labels_grouped.prop_ratio.min()

# get original proportion and counts that "determine" sampling
leading_orig_prop = remain_features_labels_grouped.loc[remain_features_labels_grouped['prop_ratio'] == min_prop_ratio_remain, 'orig_prop'].values.item()
leading_counts = remain_features_labels_grouped.loc[remain_features_labels_grouped['prop_ratio'] == min_prop_ratio_remain, 'counts'].values.item()


# get the number of instances that should be sampled from each group 
remain_features_labels_grouped['samples'] = leading_counts * remain_features_labels_grouped['orig_prop'] / leading_orig_prop 

# round down samples
remain_features_labels_grouped['samples'] = np.floor(remain_features_labels_grouped['samples']).astype(int)

# check if sample proportion is the same as original proportion
total_samples = np.sum(remain_features_labels_grouped.samples)
remain_features_labels_grouped['prop_samples'] = remain_features_labels_grouped['samples']/total_samples

# create a dictionary for 
sample_counts = remain_features_labels_grouped[['IMPACT_P_LP', 'samples']].set_index('IMPACT_P_LP')['samples'].to_dict()


# create the sampled test df 
sampled_test_set = pd.DataFrame()

for IMPACT_P_LP_value, sample_size in sample_counts.items():
    # filter the original DataFrame for the current IMPACT value
    IMPACT_P_LP_df = remain_features_labels[remain_features_labels['IMPACT_P_LP'] == IMPACT_P_LP_value]
    
    # sample the required number of rows
    new_samples = IMPACT_P_LP_df.sample(n=sample_size, replace=False, random_state=1)
    
    # append the sampled data to the final DataFrame
    sampled_test_set = pd.concat([sampled_test_set, new_samples])

# reset the index of the final sampled DataFrame
sampled_test_set.reset_index(drop=True, inplace=True)

# shuffle df
sampled_test_set = sampled_test_set.sample(replace=False, frac=1, random_state=1).reset_index(drop=True)


# save features and labels
feat_test_orig_prop = sampled_test_set.drop(columns=['P_LP', 'IMPACT_P_LP'])
labels_test_orig_prop = sampled_test_set[['var_ID', 'P_LP']]

feat_test_orig_prop.to_csv(f"variant_score/training/train_test_data/feat_test_orig_prop_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')
labels_test_orig_prop.to_csv(f"variant_score/training/train_test_data/labels_test_orig_prop_{config_vars['data_prep_date_vs']}.csv.gz", index=False, compression='gzip')




# In[13]:


len(feat_test_orig_prop)


# In[12]:


import matplotlib.pyplot as plt

plt.hist(sampled_test_set.IMPACT_num) 


# In[66]:


plt.hist(features_labels.IMPACT_num) 


# In[64]:


plt.hist(sampled_test_set.P_LP) 


# In[67]:


plt.hist(features_labels.P_LP) 


# In[33]:


plt.hist(feat_train_sub.IMPACT_num)  


# In[52]:


# find a numerical solution that keeps both the proportion of IMPACT_num equal (0.25 per group) and the 
# proportion of P_LP to 0.5


def calculate_group_no(x, p1, p2, p3, b1_max, b2_max, b3_max, b4_max, p4_max):
    """
    x = no. total instances
    p1 = no. pathogenic instances with IMPACT = 1
    p2 = no. pathogenic instances with IMPACT = 2
    p3 = no. pathogenic instances with IMPACT = 3
    """
    # b = benign, p = pathogenic
    b1 = x/4 - p1
    b2 = x/4 - p2
    b3 = x/4 - p3
    
    p4 = (0.5 * x - p1 - p2 - p3)
    b4 = x/4 - p4
    
    allowed = True
    
    # check if any number is negative
    if b1 <= 0 or b2 <= 0 or b3 <= 0 or b4 <= 0 or p4 <= 0:
        allowed = False
    
        # check if any number is above is maximum allowed value
    else: 
        if b1 > b1_max or b2 > b2_max  or b3 > b3_max  or b4 > b4_max  or p4 > p4_max:
            allowed = False
    

    return pd.DataFrame({'x': [x], 'b1': [b1], 'b2': [b2], 'b3': [b3], 'b4': [b4], 'p1': [p1], 'p2': [p2], 'p3': [p3],
                         'p4': [p4], 'allowed': [allowed]})


x_max = features_labels_grouped.size().min() * 4
b1_max = np.min([x_max/4, 51479])
p1_max = np.min([x_max/4, 357])
b2_max = np.min([x_max/4, 112283])
p2_max = np.min([x_max/4, 853])
b3_max = np.min([x_max/4, 17008])
p3_max = np.min([x_max/4, 12633])
b4_max = np.min([x_max/4, 1087])
p4_max = np.min([x_max/4, 55904])


df_sum = pd.DataFrame()


for p1 in range(0, int(p1_max), 10):
    for p2 in range(0, int(p2_max), 10):
        for p3 in range(8000, int(p3_max), 1000):
            for x in [50000, 60000, 70000, 80000, 90000]:
                
                df_row = calculate_group_no(x=x, 
                                            p1=p1,
                                            p2=p2, 
                                            p3=p3, 
                                            b1_max=b1_max, 
                                            b2_max=b2_max, 
                                            b3_max=b3_max, 
                                            b4_max=b4_max, 
                                            p4_max=p4_max)
                
                df_sum = pd.concat([df_sum, df_row])
                
                


# In[57]:


df_sum.query("allowed == True")


# In[60]:


df_sum.query("allowed == True").b3.unique()

