#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
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
import json

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

from training.helper_functions_ML import *


# In[96]:


def get_genes_for_prediction(config_dic, gene_set):
    """
    Function to return filled and scaled gene set.
    'gene_set' must be 'train', 'test', 'train_test', 'all' or 'unknown'.
    
    """
    
    ## prepare full gene set for prediction
    # load full gene set
    raw_feat = pd.read_csv(f"features/results/gene_features_{config_vars['creation_date']}.csv.gz")

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

    ## select genes for predictions    
    # machine learning test set
    if gene_set == 'test':
        feat_test = pd.read_csv(f"training/train_test_data/feat_test_{config_vars['data_prep_date']}.csv.gz")
        hgnc_ids_for_prediction = list(feat_test['hgnc_id_int'])
    
    # machine learning training set
    elif gene_set == 'train':
        feat_train = pd.read_csv(f"training/train_test_data/feat_train_{config_vars['data_prep_date']}.csv.gz")
        hgnc_ids_for_prediction = list(feat_train['hgnc_id_int'])
    
    # machine learning training and test set
    elif gene_set == 'train_test':
        feat_test = pd.read_csv(f"training/train_test_data/feat_test_{config_vars['data_prep_date']}.csv.gz")
        feat_train = pd.read_csv(f"training/train_test_data/feat_train_{config_vars['data_prep_date']}.csv.gz")
        hgnc_ids_for_prediction = list(feat_test['hgnc_id_int']) + list(feat_train['hgnc_id_int'])
    
    # all genes 
    elif gene_set == 'all':
        hgnc_ids_for_prediction = list(all_genes_scaled['hgnc_id_int'])
    
    # genes without known labels (all but training and test genes)
    elif gene_set == 'unknown':
        feat_test = pd.read_csv(f"training/train_test_data/feat_test_{config_vars['data_prep_date']}.csv.gz")
        feat_train = pd.read_csv(f"training/train_test_data/feat_train_{config_vars['data_prep_date']}.csv.gz")
        hgnc_ids_for_prediction = list(set(all_genes_scaled['hgnc_id_int']) - set(feat_test['hgnc_id_int']) - set(feat_train['hgnc_id_int']))
        
    # error in case of undefined/invalid gene set     
    else:
        raise ValueError("'gene_set' must be 'train', 'test', 'train_test', 'all' or 'unknown'.")
    
    # filter genes for prediction
    genes_for_predictions = all_genes_scaled.query("hgnc_id_int in @hgnc_ids_for_prediction")

    return genes_for_predictions


# In[56]:


def get_symbol_from_hgnc_id_int(hgnc_id_int_list):
    """
    Function that annotates a gene list of HGNC IDs with their symbol based 
    on the HGNC table of kidney-genetics on github.
    """
    # HGNC table gitHub raw URL 
    url = 'https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/B_AnnotationHGNC/results/non_alt_loci_set_coordinates.2023-11-21.csv.gz'

    # read the .csv file into a DataFrame
    hgnc_annotated = pd.read_csv(url, compression='gzip')

    # add a new column without the "HGNC:" prefix
    hgnc_annotated['hgnc_id_int'] = hgnc_annotated['hgnc_id'].str.replace('HGNC:', '')

    # convert the 'hgnc_id_without_prefix' column to integers
    hgnc_annotated['hgnc_id_int'] = pd.to_numeric(hgnc_annotated['hgnc_id_int'], downcast='integer')
    
    # create df from gene list
    genes_df = pd.DataFrame({'hgnc_id_int' : hgnc_id_int_list})
    
    # annotate with symbol
    genes_df = genes_df.merge(hgnc_annotated[['hgnc_id_int', 'symbol']], how='left', left_on='hgnc_id_int', right_on='hgnc_id_int')

    return genes_df['symbol'].tolist()


# In[42]:


def get_gene_set_from_hgnc_id_int(hgnc_id_int_list):
    """
    Function that annotates a gene list of HGNC IDs with their gene set (train, test, None).
    """
    test_ids = pd.read_csv(f"training/train_test_data/feat_test_{config_vars['data_prep_date']}.csv.gz")['hgnc_id_int'].tolist()
    train_ids = pd.read_csv(f"training/train_test_data/feat_train_{config_vars['data_prep_date']}.csv.gz")['hgnc_id_int'].tolist()
    
    def classify_gene_set(hgnc_id_int):
        if hgnc_id_int in train_ids:
            return 'train'
        elif hgnc_id_int in test_ids:
            return 'test'
        else:
            return None
    
    # create df from gene list
    genes_df = pd.DataFrame({'hgnc_id_int' : hgnc_id_int_list})
    
    # create the 'gene_set' column based on the classification function
    genes_df['gene_set'] = genes_df['hgnc_id_int'].apply(lambda x: classify_gene_set(x))
    
    return genes_df['gene_set'].tolist()


# In[97]:


def get_evidence_count_from_hgnc_id_int(hgnc_id_int_list):
    """
    Function that annotates a gene list of HGNC IDs with their evidence count.
    """
    # get positive and dispensible/negative genes
    pos_genes = pd.read_csv(f"labels/results/positive_genes_{config_vars['creation_date']}.csv.gz")
    neg_genes = pd.read_csv(f"labels/results/dispensible_genes_{config_vars['creation_date']}.csv.gz")
    
    # set evidence count of dispensible genes to -1
    neg_genes['evidence_count'] = -1

    # rowbind both
    annotated_genes = pd.concat([pos_genes[['hgnc_id_int', 'evidence_count']], neg_genes], ignore_index=True)
    
    # create df from gene list
    genes_df = pd.DataFrame({'hgnc_id_int' : hgnc_id_int_list})
    
    # annotate genes df with evidence count
    genes_df = genes_df.merge(annotated_genes[['hgnc_id_int', 'evidence_count']], how='left', left_on='hgnc_id_int', right_on='hgnc_id_int')
    
    return genes_df['evidence_count'].tolist()


# In[98]:


def make_predictions(ID, gene_set, save):
    # get configuration dictionary and results dictionary of the chosen experiment
    config_dic, results_dic = get_config_results_dics(ID=ID) 
    
    # get best parameters
    best_params = results_dic['best_params']
    
    # create classifier with best parameters    
    clf = config_dic['clf']
    
    # set estimator and best parameters
    clf.set_params(estimator=config_dic['estimator'])
    clf.set_params(**best_params)
    
    # fit classifier with training data
    clf.fit(config_dic['X_train'], config_dic['y_train'])
    
    # get gene set for prediction
    genes_for_prediction = get_genes_for_prediction(config_dic=config_dic, gene_set=gene_set)

    # get gene features as numpy arrays
    X = genes_for_prediction.drop(columns=['hgnc_id_int']).values
    X_hgnc_id_int = genes_for_prediction['hgnc_id_int']

    
    ## probability predicition
    # predict probabilities for selected genes (=> 2 dim array, probabilities sum up to 1)
    probabilities = clf.predict_proba(X)

    # get disease gene probabilities
    disease_gene_prob = probabilities[:, 1]

    # create the dataframe with Nephro Gene Score (NGS)
    NGS = pd.DataFrame({'hgnc_id_int': X_hgnc_id_int, 'NGS': disease_gene_prob})
    
    # annotate with symbols
    NGS['symbol'] = get_symbol_from_hgnc_id_int(NGS['hgnc_id_int'].tolist())
    
    # annotate with gene set
    NGS['gene_set'] = get_gene_set_from_hgnc_id_int(NGS['hgnc_id_int'].tolist())
    
    # annotate with evidence count 
    NGS['evidence_count'] = get_evidence_count_from_hgnc_id_int(NGS['hgnc_id_int'].tolist())
    
    # save csv
    if save:
        current_date = datetime.now().strftime("%Y-%m-%d")
        NGS.to_csv(f"predictions/results/NGS_predictions_ID{ID}_{gene_set}_{current_date}.csv.gz", index=False, compression='gzip')
        
    return NGS


# In[46]:


# write NCS to json file per gene

# convert to JSON with column names as keys and values in a list

# for row in np.arange(NCS_df.shape[0]):
#     json_data = NCS_df.iloc[row].to_dict()
#     hgnc_id_int = json_data['hgnc_id_int']
#     symbol = json_data['symbol']

#     # save to a file
#     with open(f'../predictions/results/json/hgnc/{hgnc_id_int}.json', 'w') as json_file:
#         json.dump(json_data, json_file)
        
#     # save to a file
#     with open(f'../predictions/results/json/symbols/{symbol}.json', 'w') as json_file:
#         json.dump(json_data, json_file)



# In[ ]:


# CONTINUE HERE: 2024-02-28
# TODO: write violin plots and boxplots as functions. write a function for creating the json files


# In[92]:


NGS = make_predictions(ID=90, gene_set='all', save=True) # TODO: add as variable in function!!!
NGS 
gene_set = 'all' # TODO: add as variable in function!!!
ID = 90 # TODO: add as variable in function!!!

## create NCS boxplots based on evidence counts
evidence_counts = [-1, 0, 1, 2, 3, 4, 5]

plt.figure(figsize=(7, 6))
for count in evidence_counts:
    subset_df = NGS[NGS['evidence_count'] == count]
    plt.boxplot(subset_df['NGS'], positions=[count], widths=0.6, showfliers=False, vert=True)

    # Add count above the boxplot
    num_genes = len(subset_df)
    if num_genes > 0:
        plt.text(count, -0.05, f'{num_genes}', ha='center', va='bottom')
        
# add the boxplot for evidence_count = nan
subset_nan = NGS[NGS['evidence_count'].isna()]
if not subset_nan.empty:
    boxplot_nan = plt.boxplot(subset_nan['NGS'], positions=[6], widths=0.6, showfliers=False, vert=True, boxprops=dict(color='red'))
    num_genes_nan = len(subset_nan)
    plt.text(6, -0.05, f'{num_genes_nan}', ha='center', va='bottom')
        
# plot the gene number label
plt.text(7, -0.05, 'no. genes', ha='center', va='bottom')

# add a horizontal line
plt.axhline(y=-0.01, color='black', linestyle='--', xmin=0, xmax=1)

        
# Set labels and title
plt.xlabel('Evidence Count')
plt.ylabel('Nephro Gene Score')
plt.title('Boxplots of Nephro Gene Score by Evidence Count')

# Customize x-axis tick labels
custom_labels = {6: 'new \n genes', -1: 'dispensible \n genes'}  # Set 'apple' label for x-position 7
tick_labels = [custom_labels.get(pos, pos) for pos in [-1, 0, 1, 2, 3, 4, 5, 6]]
plt.xticks([-1, 0, 1, 2, 3, 4, 5, 6], tick_labels)

# save plot
current_date = datetime.now().strftime("%Y-%m-%d")
plt.savefig(f"predictions/results/boxplots_NGS_by_EC_ID{ID}_{gene_set}_{current_date}.png", dpi = 300, format='png')


# Show the plot
plt.show()


# In[64]:





# In[95]:


### VIOLIN PLOT ###
evidence_counts = [-1, 0, 1, 2, 3, 4, 5]

plt.figure(figsize=(7, 6))

for count in evidence_counts:
    subset_df = NGS[NGS['evidence_count'] == count]
    
    # plot violin plot
    violin_parts = plt.violinplot(dataset=[subset_df['NGS']], positions=[count], widths=0.6, vert=True, 
                                  showextrema=False)

    # set facecolor for the violin plot
    for part in violin_parts['bodies']:
        part.set_facecolor('red')  
    
    # add count above the violin plot
    num_genes = len(subset_df)
    if num_genes > 0:
        plt.text(count, -0.05, f'{num_genes}', ha='center', va='bottom')
        
    # add median line in black
    median_val = np.median(subset_df['NGS'])
    plt.plot([count - 0.2, count + 0.2], [median_val, median_val], color='black', linestyle='-', linewidth=1)

    # add mean line in red
    mean_val = np.mean(subset_df['NGS'])
    plt.plot([count - 0.2, count + 0.2], [mean_val, mean_val], color='red', linestyle='-', linewidth=1)


# add label for median an mean
plt.plot([], [], color='black', linestyle='-', linewidth=1, label='median')
plt.plot([], [], color='red', linestyle='-', linewidth=1, label='mean')

  
# add the violin plot for evidence_count = nan
subset_nan = NGS[NGS['evidence_count'].isna()]
if not subset_nan.empty:
    violinplot_nan = plt.violinplot(dataset=[subset_nan['NGS']], positions=[6], widths=0.6, vert=True, showextrema=False)
    num_genes_nan = len(subset_nan)
    plt.text(6, -0.05, f'{num_genes_nan}', ha='center', va='bottom')

    
    # add median line in black for NaN
    median_nan = np.median(subset_nan['NGS'])
    plt.plot([6 - 0.2, 6 + 0.2], [median_nan, median_nan], color='black', linestyle='-', linewidth=1)

    # add mean line in red for NaN
    mean_nan = np.mean(subset_nan['NGS'])
    plt.plot([6 - 0.2, 6 + 0.2], [mean_nan, mean_nan], color='red', linestyle='-', linewidth=1)
    

# plot the gene number label
plt.text(7.5, -0.05, 'no. genes', ha='center', va='bottom')


# add a horizontal line
plt.axhline(y=-0.01, color='black', linestyle='--', xmin=0, xmax=1)

# set labels and title
plt.xlabel('Evidence Count')
plt.ylabel('Nephro Gene Score')
plt.title('Violinplots of Nephro Gene Score by Evidence Count')

# sustomize x-axis tick labels
custom_labels = {6: 'new \n genes', -1: 'dispensable \n genes'}  # Set 'apple' label for x-position 7
tick_labels = [custom_labels.get(pos, pos) for pos in [-1, 0, 1, 2, 3, 4, 5, 6]]
plt.xticks([-1, 0, 1, 2, 3, 4, 5, 6], tick_labels)

# add legend
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# save plot
current_date = datetime.now().strftime("%Y-%m-%d")
plt.savefig(f"predictions/results/violinplots_NGS_by_EC_ID{ID}_{gene_set}_{current_date}.png", dpi=300, format='png')

# show the plot
plt.show()


# In[ ]:


# NOTES - TO BE DELETED

# import preprocessing functions
# sys.path.append(f"{config_vars['ML_projectsdir']}{project_name}/gene_score/training")
# from functions.helper_functions_ML import *



# script_path = "/gene_score/predictions/" # TODO remove second time, as soon as paths are set to absolute path in helper_functions_ML.py



# example
a = get_genes_for_prediction(config_dic = config_dic, gene_set = 'unknown')
a 


# jan = get_symbol_from_hgnc_id_int([37133, 24086, 5, 100000000000])

get_evidence_count_from_hgnc_id_int([37133, 24086, 5, 100000000000, 36, 38])
        
    
# get_gene_set_from_hgnc_id_int([37133, 24086, 5, 100000000000, 36, 38])
    
    
# example
j = make_predictions(ID=90, gene_set='all', save=True)


# pos_genes = pd.read_csv(f"labels/results/positive_genes_{config_vars['creation_date']}.csv.gz")
# neg_genes = pd.read_csv(f"labels/results/dispensible_genes_{config_vars['creation_date']}.csv.gz")

# # # merge 
# # NCS_df = pd.merge(NCS_df, pos_genes[['hgnc_id_int', 'evidence_count']], how='left', left_on='hgnc_id_int', right_on='hgnc_id_int')

# # # set evidence_count to -1 for dispensible genes
# # mask = NCS_df['hgnc_id_int'].isin(neg_genes['hgnc_id_int'])
# # NCS_df.loc[mask, 'evidence_count'] = -1

# pos_genes
# neg_genes['evidence_count'] = -1
# annotated_genes = pd.concat([pos_genes[['hgnc_id_int', 'evidence_count']], neg_genes], ignore_index=True)
# annotated_genes.shape

    

    
    
# # create final model to predict values for all genes
# from xgboost import XGBClassifier

# # define training ID of the finally chosen model
# ID = 90 # TODO: change, if needed

# # get config_dic and results_dic
# config_dic, results_dic = get_config_results_dics(ID=ID)    

# # get best params
# best_params = results_dic['best_params']

# ## create XGBoost classifier with best params
# estimator = None
# clf = XGBClassifier(random_state=1, 
#                     booster='gbtree',
#                     gamma=best_params['gamma'],
#                     learning_rate=best_params['gamma'],
#                     max_depth=best_params['max_depth'],
#                     n_estimators=best_params['n_estimators'],
#                     reg_alpha=best_params['reg_alpha'],
#                     reg_lambda=best_params['reg_lambda'],
#                     subsample=best_params['subsample']
#                    )

# # fit classifier with training data
# clf.fit(config_dic['X_train'], config_dic['y_train'])


# ## prepare full gene set for prediction
# # load gene set
# raw_feat = pd.read_csv(f"../features/results/gene_features_{config_vars['creation_date']}.csv.gz")

# # select only features that were also used in the training process
# all_genes_df = raw_feat[['hgnc_id_int'] + config_dic['features']] # TODO: if drop features, add!

# # fill missing values
# all_genes_filled = fill_missing_vals(all_genes_df, config_dic['model'])

# # get features that should (not) be scaled and scaling method
# omit_scaling_features = config_dic['omit_scaling_features']
# scaling_features = [i for i in config_dic['features'] if i not in omit_scaling_features] # features that should be scaled
    
# scaling = config_dic['scaling']
    
# # scale features
# if scaling == 'standard':
#     # create StandardScaler
#     stand_scal = StandardScaler()
    
#     # scale features
#     all_genes_scaled = all_genes_filled.copy()
#     all_genes_scaled[scaling_features] = stand_scal.fit_transform(all_genes_scaled[scaling_features])

        
# elif scaling == 'robust':
#     # create RobustScaler
#     rob_scal = RobustScaler(with_centering=True, with_scaling=True)
    
#     # scale features
#     all_genes_scaled = all_genes_filled.copy()
#     all_genes_scaled[scaling_features] = rob_scal.fit_transform(all_genes_scaled[scaling_features])
        

# # get all scaled gene features as numpy arrays
# X_all = all_genes_scaled.drop(columns=['hgnc_id_int']).values
# X_all_hgnc_id_int = all_genes_scaled['hgnc_id_int']

# ## probability predicition
# # predict probabilities for all genes (=> 2 dim array, probabilities sum up to 1)
# probabilities = clf.predict_proba(X_all)

# # get disease gene probabilities
# disease_gene_prob = probabilities[:, 1]

# # create the NCS df
# NCS_df = pd.DataFrame({'hgnc_id_int': X_all_hgnc_id_int, 'NCS': disease_gene_prob})


# ## get annotated HGNC table from kidney-genetics
# # gitHub raw URL 
# url = 'https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/B_AnnotationHGNC/results/non_alt_loci_set_coordinates.2023-11-21.csv.gz'

# # read the CSV file into a DataFrame
# hgnc_annotated = pd.read_csv(url, compression='gzip')

# # add a new column without the "HGNC:" prefix
# hgnc_annotated['hgnc_id_int'] = hgnc_annotated['hgnc_id'].str.replace('HGNC:', '')

# # convert the 'hgnc_id_without_prefix' column to integers
# hgnc_annotated['hgnc_id_int'] = pd.to_numeric(hgnc_annotated['hgnc_id_int'], downcast='integer')

# # annotate NCS_df with symbol
# NCS_df = pd.merge(NCS_df, hgnc_annotated[['hgnc_id_int', 'symbol']], how='left', left_on='hgnc_id_int', right_on='hgnc_id_int')

# # write csv
# current_date = datetime.now().strftime("%Y-%m-%d")
# # NCS_df.to_csv(f"results/NCS_predictions_ID{ID}_{current_date}.csv.gz", index=False, compression='gzip')


# In[ ]:


# NOTES - TO BE DELETED


# ## create NCS boxplots based on evidence counts
# evidence_counts = [-1, 0, 1, 2, 3, 4, 5]

# plt.figure(figsize=(7, 6))
# for count in evidence_counts:
#     subset_df = NCS_df[NCS_df['evidence_count'] == count]
#     plt.boxplot(subset_df['NCS'], positions=[count], widths=0.6, showfliers=False, vert=True)

#     # Add count above the boxplot
#     num_genes = len(subset_df)
#     if num_genes > 0:
#         plt.text(count, -0.05, f'{num_genes}', ha='center', va='bottom')
        
# # add the boxplot for evidence_count = nan
# subset_nan = NCS_df[NCS_df['evidence_count'].isna()]
# if not subset_nan.empty:
#     boxplot_nan = plt.boxplot(subset_nan['NCS'], positions=[6], widths=0.6, showfliers=False, vert=True, boxprops=dict(color='red'))
#     num_genes_nan = len(subset_nan)
#     plt.text(6, -0.05, f'{num_genes_nan}', ha='center', va='bottom')
        
# # plot the gene number label
# plt.text(7, -0.05, 'no. genes', ha='center', va='bottom')

# # add a horizontal line
# plt.axhline(y=-0.01, color='black', linestyle='--', xmin=0, xmax=1)

        
# # Set labels and title
# plt.xlabel('Evidence Count')
# plt.ylabel('NCS')
# plt.title('Boxplots of NCS by Evidence Count')

# # Customize x-axis tick labels
# custom_labels = {6: 'new \n genes', -1: 'dispensible \n genes'}  # Set 'apple' label for x-position 7
# tick_labels = [custom_labels.get(pos, pos) for pos in [-1, 0, 1, 2, 3, 4, 5, 6]]
# plt.xticks([-1, 0, 1, 2, 3, 4, 5, 6], tick_labels)

# # save plot
# # plt.savefig(f"results/boxplots_NCS_by_EC_ID{ID}_{current_date}.png", dpi = 300, format='png')


# # Show the plot
# plt.show()

