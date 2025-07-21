#!/usr/bin/env python
# coding: utf-8

#### Results and final model ####
# loads the final model, analyses its results/feature importance (SHAP values)
# converted from .ipynb


# import basic modules
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import pandas as pd
import sys
import yaml

# import sklearn functions
from sklearn.tree import plot_tree
from sklearn.metrics import roc_auc_score
from sklearn.tree import DecisionTreeClassifier, export_text
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LogisticRegression


# set options
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

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

score = 'vs'
score_string, id_string, label_string = get_score_specific_args(score)


# In[3]:


config_dir = "variant_score/training/config_files"

# read latest config table
print(get_latest_ML_config_training_file(config_dir, score))
pd.read_csv(get_latest_ML_config_training_file(config_dir, score))


# In[4]:


ID = 3 # finally chosen LogReg model
config_dic, results_dic = get_config_results_dics(ID=ID, score=score)    

print(results_dic['best_params'])
results_dic['cv_results']['mean_test_score'].max()


# In[6]:


feat_train = pd.read_csv(f"variant_score/training/train_test_data/feat_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz")
labels_train = pd.read_csv(f"variant_score/training/train_test_data/labels_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz")

labels_feat_train = pd.merge(labels_train, feat_train, on='var_ID', how='inner') 
labels_feat_train[['var_ID', 'P_LP'] + config_dic['features']].head(10)


# In[9]:


# plot 2D heatmap with fixed other hyperparamters (to best values)    
plot_2D_heatmap_fixed_params(ID=ID, 
                            cv_results=results_dic['cv_results'], 
                            best_params = results_dic['best_params'],
                            param1='penalty', 
                            param2='C', 
                            figsize=(15,10), 
                            save=False, 
                            show=True)


# In[7]:


clf = results_dic['best_classifier']
clf.fit(config_dic['X_train'], config_dic['y_train'])


# In[8]:


# Extract coefficients and intercept for LogReg model
coefficients = clf.coef_
intercept = clf.intercept_

print(config_dic['features'])
print(coefficients)
print(intercept)

pd.DataFrame({'feature':config_dic['features'], 'coeff': coefficients[0] }).sort_values("coeff")


# In[ ]:





# In[9]:


feat = pd.DataFrame(config_dic['X_train'])
feat.columns = config_dic['features']
feat_labels = pd.concat([feat, pd.DataFrame({'P_LP': config_dic['y_train']})], axis=1)
feat_labels['y_pred'] = clf.predict_proba(config_dic['X_train'])[:, -1]
feat_labels['log_proba'] = clf.predict_log_proba(config_dic['X_train'])[:, -1]
feat_labels


# In[10]:


sampled_df = feat_labels.sample(n=1000)
color_map = {1: 'red', 0: 'blue'}

plt.scatter(sampled_df['log_proba'].values, sampled_df['y_pred'].values,
           c=sampled_df['P_LP'].map(color_map),
           cmap=plt.cm.bwr,
           s=1)
plt.xlabel('log_proba')
plt.ylabel('pred')
plt.xlim([-10, sampled_df['log_proba'].values.max()])


# In[23]:


# SHAP values
import shap
import warnings
from sklearn.preprocessing import StandardScaler


warnings.filterwarnings("ignore", message="is_sparse is deprecated and will be removed in a future version.")

ID = 3
config_dic, results_dic = get_config_results_dics(ID=ID, score='vs') 
features = config_dic['features']

random.seed(1)
sample_indices = random.sample(range(0, config_dic['X_train'].shape[0]), 900) # CAVE: low sample size => result very instable!!!

# X = config_dic['X_train'][sample_indices, :]
X = config_dic['X_train'] #[sample_indices, :]


# initialize StandardScaler
# scaler = StandardScaler()

# fit the scaler to your training data and transform it
# X_scaled = scaler.fit_transform(X)

# X_lab = pd.DataFrame(X_scaled, columns=features)
X_lab = pd.DataFrame(X, columns=features)

# y = config_dic['y_train'][sample_indices]
y = config_dic['y_train']
clf = results_dic['best_classifier']
clf.random_state = 1 

model = clf.fit(X_lab, y)

explainer = shap.Explainer(model.predict, X_lab)  
# explainer = shap.Explainer(model, X_lab)  # 38

# explainer = shap.TreeExplainer(model, X_lab ) 
shap_values = explainer(X_lab)

# # dump SHAP values in pickle - DONT DUMP IN PICKLE - USE JSON
# with open(f"training/feature_importance/SHAP_values/values/SHAP_values_ID{ID}_{datetime.today().strftime('%Y-%m-%d')}.pkl", 'wb') as file:
#     pickle.dump(shap_values, file)


# In[17]:


shap.summary_plot(shap_values, max_display=20, show=True, plot_size=[10, 8]) 

