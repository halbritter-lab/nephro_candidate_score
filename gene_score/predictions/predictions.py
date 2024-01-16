#!/usr/bin/env python
# coding: utf-8

# In[2]:


os.environ['CONFIG_FILE'] = '/fast/work/users/rankn_c/halbritter/nephro_candidate_score/gene_score/training/config_NCS.yml' # TODO: change


# In[ ]:





# In[3]:


# import basic modules
import numpy as np
import pandas as pd
import matplotlib as plt
import yaml
import os
import sys

# set options
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"
script_path = "/gene_score/predictions/"

# read configs
CONFIG_FILE = os.getenv('CONFIG_FILE')

with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}{script_path}")

# import preprocessing functions
sys.path.append(f"{config_vars['ML_projectsdir']}{project_name}/gene_score/training")

from functions.helper_functions_ML import *

script_path = "/gene_score/predictions/" # TODO remove second time, as soon as paths are set to absolute path in helper_functions_ML.py


# In[16]:


# import sys
# sys.path.append(f"{config_vars['ML_projectsdir']}{project_name}/gene_score/training")

# from helper_functions_ML import *


# In[4]:


# create final model to predict values for all genes
from xgboost import XGBClassifier

# define training ID of the finally chosen model
ID = 90 # TODO: change, if needed

# get config_dic and results_dic
config_dic, results_dic = get_config_results_dics(ID=ID)    

# get best params
best_params = results_dic['best_params']

## create XGBoost classifier with best params
estimator = None
clf = XGBClassifier(random_state=1, 
                    booster='gbtree',
                    gamma=best_params['gamma'],
                    learning_rate=best_params['gamma'],
                    max_depth=best_params['max_depth'],
                    n_estimators=best_params['n_estimators'],
                    reg_alpha=best_params['reg_alpha'],
                    reg_lambda=best_params['reg_lambda'],
                    subsample=best_params['subsample']
                   )

# fit classifier with training data
clf.fit(config_dic['X_train'], config_dic['y_train'])


## prepare full gene set for prediction
# load gene set
raw_feat = pd.read_csv(f"../features/results/gene_features_{config_vars['creation_date']}.csv.gz")

# select only features that were also used in the training process
all_genes_df = raw_feat[['hgnc_id_int'] + config_dic['features']] # TODO: if drop features, add!

# fill missing values
all_genes_filled = fill_missing_vals(all_genes_df, config_dic['model'])

# get features that should (not) be scaled and scaling method
omit_scaling_features = config_dic['omit_scaling_features']
scaling_features = [i for i in config_dic['features'] if i not in omit_scaling_features] # features that should be scaled
    
scaling = config_dic['scaling']
    
# scale features
if scaling == 'standard':
    # create StandardScaler
    stand_scal = StandardScaler()
    
    # scale features
    all_genes_scaled = all_genes_filled.copy()
    all_genes_scaled[scaling_features] = stand_scal.fit_transform(all_genes_scaled[scaling_features])

        
elif scaling == 'robust':
    # create RobustScaler
    rob_scal = RobustScaler(with_centering=True, with_scaling=True)
    
    # scale features
    all_genes_scaled = all_genes_filled.copy()
    all_genes_scaled[scaling_features] = rob_scal.fit_transform(all_genes_scaled[scaling_features])
        

# get all scaled gene features as numpy arrays
X_all = all_genes_scaled.drop(columns=['hgnc_id_int']).values
X_all_hgnc_id_int = all_genes_scaled['hgnc_id_int']

## probability predicition
# predict probabilities for all genes (=> 2 dim array, probabilities sum up to 1)
probabilities = clf.predict_proba(X_all)

# get disease gene probabilities
disease_gene_prob = probabilities[:, 1]

# create the NCS df
NCS_df = pd.DataFrame({'hgnc_id_int': X_all_hgnc_id_int, 'NCS': disease_gene_prob})


## get annotated HGNC table from kidney-genetics
# gitHub raw URL 
url = 'https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/B_AnnotationHGNC/results/non_alt_loci_set_coordinates.2023-11-21.csv.gz'

# read the CSV file into a DataFrame
hgnc_annotated = pd.read_csv(url, compression='gzip')

# add a new column without the "HGNC:" prefix
hgnc_annotated['hgnc_id_int'] = hgnc_annotated['hgnc_id'].str.replace('HGNC:', '')

# convert the 'hgnc_id_without_prefix' column to integers
hgnc_annotated['hgnc_id_int'] = pd.to_numeric(hgnc_annotated['hgnc_id_int'], downcast='integer')

# annotate NCS_df with symbol
NCS_df = pd.merge(NCS_df, hgnc_annotated[['hgnc_id_int', 'symbol']], how='left', left_on='hgnc_id_int', right_on='hgnc_id_int')

# write csv
current_date = datetime.now().strftime("%Y-%m-%d")
NCS_df.to_csv(f"results/NCS_df_ID{ID}_{current_date}.csv.gz", index=False, compression='gzip')


# In[5]:


NCS_df


# In[6]:


## merge with evidence count
# load positive and dispensible genes
pos_genes = pd.read_csv(f"../labels/results/positive_genes_{config_vars['creation_date']}.csv.gz")
neg_genes = pd.read_csv(f"../labels/results/dispensible_genes_{config_vars['creation_date']}.csv.gz")

# merge 
NCS_df = pd.merge(NCS_df, pos_genes[['hgnc_id_int', 'evidence_count']], how='left', left_on='hgnc_id_int', right_on='hgnc_id_int')

# set evidence_count to -1 for dispensible genes
mask = NCS_df['hgnc_id_int'].isin(neg_genes['hgnc_id_int'])
NCS_df.loc[mask, 'evidence_count'] = -1


# In[7]:


## create NCS boxplots based on evidence counts
evidence_counts = [-1, 0, 1, 2, 3, 4, 5]

plt.figure(figsize=(8, 7))
for count in evidence_counts:
    subset_df = NCS_df[NCS_df['evidence_count'] == count]
    plt.boxplot(subset_df['NCS'], positions=[count], widths=0.6, showfliers=False, vert=True)

    # Add count above the boxplot
    num_genes = len(subset_df)
    if num_genes > 0:
        plt.text(count, -0.045, f'{num_genes}', ha='center', va='bottom')
        
# add the boxplot for evidence_count = nan
subset_nan = NCS_df[NCS_df['evidence_count'].isna()]
if not subset_nan.empty:
    boxplot_nan = plt.boxplot(subset_nan['NCS'], positions=[6], widths=0.6, showfliers=False, vert=True, boxprops=dict(color='red'))
    num_genes_nan = len(subset_nan)
    plt.text(6, -0.045, f'{num_genes_nan}', ha='center', va='bottom')
        
# plot the gene number label
plt.text(7, -0.045, 'no. genes', ha='center', va='bottom')

# add a horizontal line
plt.axhline(y=-0.01, color='black', linestyle='--', xmin=0, xmax=1)

        
# Set labels and title
plt.xlabel('Evidence Count')
plt.ylabel('NCS')
plt.title('Boxplots of NCS by Evidence Count')

# Customize x-axis tick labels
custom_labels = {6: 'new \n genes', -1: 'dispensible \n genes'}  # Set 'apple' label for x-position 7
tick_labels = [custom_labels.get(pos, pos) for pos in [-1, 0, 1, 2, 3, 4, 5, 6]]
plt.xticks([-1, 0, 1, 2, 3, 4, 5, 6], tick_labels)

# save plot
plt.savefig(f"results/boxplots_NCS_by_EC_ID{ID}_{current_date}.png", dpi = 300, format='png')


# Show the plot
plt.show()


# In[ ]:




