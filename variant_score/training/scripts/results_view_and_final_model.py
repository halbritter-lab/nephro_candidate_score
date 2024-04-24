#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
os.environ['CONFIG_FILE'] = '/fast/work/users/rankn_c/halbritter/nephro_candidate_score/gene_score/training/config_NCS.yml' # TODO: change
# TODO: change


# In[3]:


# import basic modules
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import yaml

# import sklearn functions
from sklearn.tree import plot_tree
from sklearn.metrics import roc_auc_score
from sklearn.tree import DecisionTreeClassifier, export_text

# set options
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

# get config file
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

# import preprocessing functions
from helper_functions_ML_vs import *


# In[28]:


config_dir = "training/config_files"

# read latest config table
print(get_latest_ML_config_training_file(config_dir))
pd.read_csv(get_latest_ML_config_training_file(config_dir))


# In[29]:


ID = 1
config_dic, results_dic = get_config_results_dics(ID=ID)    

print(results_dic['best_params'])
results_dic['cv_results']['mean_test_score'].max()


# In[36]:


# plot 2D heatmap with fixed other hyperparamters (to best values)    
ID = 1
# config_dic, results_dic = get_config_results_dics(ID=ID)  


plot_2D_heatmap_fixed_params(ID=ID, 
                            cv_results=results_dic['cv_results'], 
                            best_params = results_dic['best_params'],
                            param1='max_depth', 
                            param2='min_samples_split', 
                            figsize=(15,10), 
                            save=False, 
                            show=True)


# In[46]:


results_dic['cv_results']


# In[33]:


# Plot decision tree
X = config_dic['X_train']
y = config_dic['y_train']

# Train a decision tree classifier
clf = results_dic['best_classifier']
clf.fit(X, y)

# Plot the decision tree
plt.figure(figsize=(10,8))
plot_tree(clf, max_depth=5, filled=True, feature_names=config_dic['features'], class_names=['(Likely) Benign','(Likely) Patho'],
         fontsize=7)
plt.show()


# In[ ]:





# In[32]:


# feature importance of best classifier
pd.DataFrame({'feature': config_dic['features'], 'feat_imp' : clf.feature_importances_}).sort_values(by='feat_imp', ascending=False)


# In[44]:


##  train DT on full training set
# define maximum depth
# TODO: less deep or deeper?
max_depth = 5 # instead of 9, as less deep DT perform already very well, use less deep DT

# define param dic
params = {'criterion': 'entropy',
          'max_depth': max_depth, 
          'min_samples_leaf': 10,
          'min_samples_split': 150}

# create classifier with best parameters    
clf = DecisionTreeClassifier(random_state=1)    

# set parameters
clf.set_params(**params)

# fit classifier
clf.fit(config_dic['X_train'], config_dic['y_train'])

# predict training set (=> 2 dim array, probabilities sum up to 1)
probabilities = clf.predict_proba(config_dic['X_train'])

# get training AUC
training_auc = roc_auc_score(config_dic['y_train'], probabilities[:, 1])  # Use the probabilities of the positive class
print("Training AUC:", training_auc)



# In[55]:


# Plot the decision tree
plt.figure(figsize=(14,10))
plot_tree(clf, filled=True, feature_names=config_dic['features'], class_names=['(Likely) Benign','(Likely) Patho'],
         fontsize=7)
plt.show()


# In[51]:


# export the decision tree to text format
tree_rules = export_text(clf, feature_names=config_dic['features'])
print(tree_rules)


# In[53]:


# feature importance
pd.DataFrame({'feature': config_dic['features'], 'feat_imp' : clf.feature_importances_}).sort_values(by='feat_imp', ascending=False)


# In[ ]:




