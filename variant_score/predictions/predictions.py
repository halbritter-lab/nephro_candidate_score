#!/usr/bin/env python
# coding: utf-8

#### Predictions for N-VS model ####
# determines AUC for different IMPACT groups and overall AUC
# converted from .ipynb


# import basic modules
import numpy as np
import pandas as pd
import yaml
import os
# import vcf
import gzip
import re
from datetime import datetime
import matplotlib.pyplot as plt
import sys
from sklearn.metrics import roc_auc_score


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

# append path where common functions are located
sys.path.append(f"{config_vars['ML_projectsdir']}{project_name}")

# import common functions
from common_functions.training_helper_functions import *

score = 'vs'
score_string, id_string, label_string = get_score_specific_args(score)


# In[5]:


import warnings

# Suppress the specific FutureWarning
warnings.filterwarnings('ignore', category=FutureWarning, message='is_sparse is deprecated and will be removed in a future version')


# In[6]:


def get_variants_for_prediction(config_dic, variant_set, IMPACT_prop_training, test_prop, IMPACT_group, score):
    """
    Function to return filled and scaled variant set.
    'variant_set' must be 'train', 'test', 'train_test'.
    
    """
    
    score_string, id_string, label_string = get_score_specific_args(score)

    # TODO: add option for IMPACT_prop_training == False
    
    # load all features and labels
    raw_feat_labels = pd.read_csv(f"{score_string}/features_labels/results/features_labels_{config_vars['data_prep_date_vs']}.csv.gz") #TODO: customize '_vs'

    # select only features that were also used in the training process
    all_variants_df = raw_feat_labels[['var_ID'] + config_dic['features']] # TODO: if drop features, add!

    # fill missing values
    all_variants_filled = fill_missing_vals(all_variants_df, config_dic['model'], score)
    
    # get training set to calculate mean and std for scaling    
    if IMPACT_prop_training:     # TODO: add option for IMPACT_prop_training == False
        feat_train = pd.read_csv(f"{score_string}/training/train_test_data/feat_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz") #TODO: customize '_vs'
        
    var_IDs_train = list(feat_train['var_ID'])
    
    # get features that should (not) be scaled and scaling method
    omit_scaling = get_features_from_groups(config_dic['omit_scaling_features'], feat_train, score)
    scaling_features = [i for i in config_dic['features'] if i not in omit_scaling]
    
    scaling = config_dic['scaling']
    print("Scaling:", scaling)

    if scaling == 'standard':
        scaler = StandardScaler()
    if scaling == 'robust':
        scaler = RobustScaler(with_centering=True, with_scaling=True)
    
    # fit the scaler only on the training data to compute mean and std
    scaler.fit(all_variants_filled.query("var_ID in @var_IDs_train")[scaling_features])
               
    # transform all variants (with mean and std from training set)
    all_variants_scaled = all_variants_filled.copy()
    all_variants_scaled[scaling_features] = scaler.transform(all_variants_scaled[scaling_features])

    ## select variants for predictions    #TODO: change!!
    # machine learning test set
    if variant_set == 'test':
        
        if test_prop == "balanced":
            feat_test = pd.read_csv(f"{score_string}/training/train_test_data/feat_test_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz") #TODO: customize '_vs'
        elif test_prop == "original":
            feat_test = pd.read_csv(f"{score_string}/training/train_test_data/feat_test_orig_prop_{config_vars['data_prep_date_vs']}.csv.gz") #TODO: customize '_vs'
        elif test_prop == "all_remaining":
            feat_test = pd.read_csv(f"{score_string}/training/train_test_data/feat_test_remain_IMPACT_{IMPACT_group}_{config_vars['data_prep_date_vs']}.csv.gz") #TODO: customize '_vs'
        else:       
            raise ValueError("'test_prop' must be TODO.")

        var_IDs_for_prediction = list(feat_test['var_ID'])
        
        
    # machine learning training set
    elif variant_set == 'train':
        feat_train = pd.read_csv(f"{score_string}/training/train_test_data/feat_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz") #TODO: customize '_vs'
        var_IDs_for_prediction = list(feat_train['var_ID'])
    
    # machine learning training and test set
    elif variant_set == 'train_test':
        feat_test = pd.read_csv(f"{score_string}/training/train_test_data/feat_test_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz") #TODO: customize '_vs'
        feat_train = pd.read_csv(f"{score_string}/training/train_test_data/feat_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz") #TODO: customize '_vs'
        var_IDs_for_prediction = list(feat_test['var_ID']) + list(feat_train['var_ID'])
    
    # error in case of undefined/invalid variant set     
    else:
        raise ValueError("'variant_set' must be 'train', 'test' or 'train_test'.")
    
    # filter variants for prediction
    variants_for_predictions = all_variants_scaled.query("var_ID in @var_IDs_for_prediction")

    return variants_for_predictions



# In[18]:


def make_predictions(ID, variant_set, IMPACT_prop_training, test_prop, IMPACT_group, save):
    # get configuration dictionary and results dictionary of the chosen experiment
    config_dic, results_dic = get_config_results_dics(ID=ID, score='vs') 
    
    # get best parameters
    best_params = results_dic['best_params']
    
    # create classifier with best parameters    
    clf = config_dic['clf']
    
    # set estimator and best parameters
    if config_dic['estimator']:
        clf.set_params(estimator=config_dic['estimator'])
    clf.set_params(**best_params)
    
    # fit classifier with training data
    clf.fit(config_dic['X_train'], config_dic['y_train'])

    # Print the coefficients
    print("Coefficients:", clf.coef_)
    
    # Print the intercept
    print("Intercept:", clf.intercept_)

    
    # get variants set for prediction
    variants_for_prediction = get_variants_for_prediction(config_dic=config_dic, 
                                                          variant_set=variant_set,
                                                          IMPACT_prop_training=IMPACT_prop_training,
                                                          test_prop=test_prop, 
                                                          IMPACT_group=IMPACT_group,
                                                          score='vs'
                                                         )
    print(len(variants_for_prediction))

    # get variants features as numpy arrays
    X = variants_for_prediction.drop(columns=['var_ID']).values
    X_var_ID = variants_for_prediction['var_ID']
    
    
    ## probability predicition
    # predict probabilities for selected variants (=> 2 dim array, probabilities sum up to 1)
    probabilities = clf.predict_proba(X)

    # get probabilities for variants beeing (likely) pathogenic
    P_LP_prob = probabilities[:, 1]

    # create the dataframe with Nephro Gene Score (NGS)
    NVS = pd.DataFrame({'var_ID': X_var_ID, 'NVS': P_LP_prob})
    
    # save csv
    if save:
        current_date = datetime.now().strftime("%Y-%m-%d")
        NVS.to_csv(f"variant_score/predictions/results/NVS_predictions_ID{ID}_{gene_set}_{current_date}.csv.gz", index=False, compression='gzip')
        
    return NVS


# In[1]:


features_labels = pd.read_csv(f"variant_score/features_labels/results/features_labels_{config_vars['data_prep_date_vs']}.csv.gz")



# In[20]:


NVS_train = make_predictions(ID=ID, variant_set='train', 
                             IMPACT_prop_training=True,  
                             test_prop=None, 
                             IMPACT_group=None,
                             save=False)
# calculate AUC
# labels_train_sub = pd.read_csv("training/train_test_data/labels_train_IMPACT_prop_2024-05-15.csv.gz")
NVS_train = pd.merge(NVS_train, features_labels[['var_ID', 'P_LP']], how='left', left_on='var_ID', right_on='var_ID')
roc_auc_score(NVS_train['P_LP'], NVS_train['NVS'])


# In[116]:


NVS_test_1 = make_predictions(ID=ID, 
                              variant_set='test', 
                              IMPACT_prop_training=True,  
                              test_prop='all_remaining', 
                              IMPACT_group=1,
                              save=False)

NVS_test_1 = pd.merge(NVS_test_1, features_labels[['var_ID', 'P_LP']], how='left', left_on='var_ID', right_on='var_ID')

roc_auc_score(NVS_test_1['P_LP'], NVS_test_1['NVS'])


# In[107]:


NVS_test_2 = make_predictions(ID=ID, 
                              variant_set='test', 
                              IMPACT_prop_training=True,  
                              test_prop='all_remaining', 
                              IMPACT_group=2,
                              save=False)

NVS_test_2 = pd.merge(NVS_test_2, features_labels[['var_ID', 'P_LP']], how='left', left_on='var_ID', right_on='var_ID')

roc_auc_score(NVS_test_2['P_LP'], NVS_test_2['NVS'])


# In[110]:


NVS_test_3 = make_predictions(ID=ID, 
                              variant_set='test', 
                              IMPACT_prop_training=True,  
                              test_prop='all_remaining', 
                              IMPACT_group=3,
                              save=False)

NVS_test_3 = pd.merge(NVS_test_3, features_labels[['var_ID', 'P_LP']], how='left', left_on='var_ID', right_on='var_ID')

roc_auc_score(NVS_test_3['P_LP'], NVS_test_3['NVS'])


# In[112]:


NVS_test_4 = make_predictions(ID=ID, 
                              variant_set='test', 
                              IMPACT_prop_training=True,  
                              test_prop='all_remaining', 
                              IMPACT_group=4,
                              save=False)

NVS_test_4 = pd.merge(NVS_test_4, features_labels[['var_ID', 'P_LP']], how='left', left_on='var_ID', right_on='var_ID')

roc_auc_score(NVS_test_4['P_LP'], NVS_test_4['NVS'])


# In[113]:


NVS_test_balanced = make_predictions(ID=ID, 
                              variant_set='test', 
                              IMPACT_prop_training=True,  
                              test_prop='balanced', 
                              IMPACT_group=None,
                              save=False)

NVS_test_balanced = pd.merge(NVS_test_balanced, features_labels[['var_ID', 'P_LP']], how='left', left_on='var_ID', right_on='var_ID')

roc_auc_score(NVS_test_balanced['P_LP'], NVS_test_balanced['NVS'])


# In[114]:


NVS_test_orig_prop = make_predictions(ID=ID, 
                              variant_set='test', 
                              IMPACT_prop_training=True,  
                              test_prop='original', 
                              IMPACT_group=None,
                              save=False)

NVS_test_orig_prop = pd.merge(NVS_test_orig_prop, features_labels[['var_ID', 'P_LP']], how='left', left_on='var_ID', right_on='var_ID')

roc_auc_score(NVS_test_orig_prop['P_LP'], NVS_test_orig_prop['NVS'])


# In[122]:


NVS_test_orig_prop = make_predictions(ID=ID, 
                              variant_set='test', 
                              IMPACT_prop_training=True,  
                              test_prop='original', 
                              IMPACT_group=None,
                              save=False)

NVS_test_orig_prop = pd.merge(NVS_test_orig_prop, features_labels[['var_ID', 'P_LP', 'IMPACT_num']], how='left', left_on='var_ID', right_on='var_ID')

NVS_test_orig_prop1 = NVS_test_orig_prop.query("IMPACT_num == 1")
print(roc_auc_score(NVS_test_orig_prop1['P_LP'], NVS_test_orig_prop1['NVS']))
print(len(NVS_test_orig_prop1))

NVS_test_orig_prop2 = NVS_test_orig_prop.query("IMPACT_num == 2")
print(roc_auc_score(NVS_test_orig_prop2['P_LP'], NVS_test_orig_prop2['NVS']))
print(len(NVS_test_orig_prop2))


NVS_test_orig_prop3 = NVS_test_orig_prop.query("IMPACT_num == 3")
print(roc_auc_score(NVS_test_orig_prop3['P_LP'], NVS_test_orig_prop3['NVS']))
print(len(NVS_test_orig_prop3))


NVS_test_orig_prop4 = NVS_test_orig_prop.query("IMPACT_num == 4")
print(roc_auc_score(NVS_test_orig_prop4['P_LP'], NVS_test_orig_prop4['NVS']))
print(len(NVS_test_orig_prop4))

all_original = pd.DataFrame()
all_original = pd.concat([all_original, NVS_test_orig_prop1])
all_original = pd.concat([all_original, NVS_test_orig_prop2])
all_original = pd.concat([all_original, NVS_test_orig_prop3])
all_original = pd.concat([all_original, NVS_test_orig_prop4])

print(len(all_original))

print(roc_auc_score(all_original['P_LP'], all_original['NVS']))


# In[ ]:


# Notes: https://discourse.datamethods.org/t/auc-of-each-subgroup-is-smaller-than-overall-auc/3264


# In[6]:


pd.read_csv(f"{score_string}/training/train_test_data/feat_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz") 
        

