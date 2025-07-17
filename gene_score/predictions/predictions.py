#!/usr/bin/env python
# coding: utf-8

#### N-GS predictions and results analysis ####
# converted from .ipynb

# TODO: get_genes_for_prediction can currently not be used for PCA reduced predictions


# import basic modules
import numpy as np
import pandas as pd
import matplotlib as plt
import yaml
import os
import sys
import json
from sklearn.metrics import roc_auc_score


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
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}")

#from training.helper_functions_ML import *


# In[6]:


# import common functions
from common_functions.training_helper_functions import *


# In[ ]:





# In[17]:


def get_genes_for_prediction(config_dic, gene_set):
    """
    Function to return filled and scaled gene set.
    'gene_set' must be 'train', 'test', 'train_test', 'all' or 'unknown'.
    
    """
    
    ## prepare full gene set for prediction
    # load full gene set
    raw_feat = pd.read_csv(f"gene_score/features/results/gene_features_{config_vars['creation_date_gs']}.csv.gz")

    # select only features that were also used in the training process
    all_genes_df = raw_feat[['hgnc_id_int'] + config_dic['features']] # TODO: if drop features, add!

    # fill missing values
    all_genes_filled = fill_missing_vals(all_genes_df, config_dic['model'], score='gs')
    
    # get training set to calculate mean and std for scaling
    feat_train = pd.read_csv(f"gene_score/training/train_test_data/feat_train_{config_vars['data_prep_date_gs']}.csv.gz")
    hgnc_ids_train = list(feat_train['hgnc_id_int'])
    
    # get features that should (not) be scaled and scaling method
    omit_scaling = get_features_from_groups(config_dic['omit_scaling_features'], feat_train, score='gs')
    scaling_features = [i for i in config_dic['features'] if i not in omit_scaling]
    
    scaling = config_dic['scaling']

    if scaling == 'standard':
        scaler = StandardScaler()
    if scaling == 'robust':
        scaler = RobustScaler(with_centering=True, with_scaling=True)
    
    # fit the scaler only on the training data to compute mean and std
    scaler.fit(all_genes_filled.query("hgnc_id_int in @hgnc_ids_train")[scaling_features])
               
    # transform all genes (with mean and std from training set)
    all_genes_scaled = all_genes_filled.copy()
    all_genes_scaled[scaling_features] = scaler.transform(all_genes_scaled[scaling_features])

    ## select genes for predictions    
    # machine learning test set
    if gene_set == 'test':
        feat_test = pd.read_csv(f"gene_score/training/train_test_data/feat_test_{config_vars['data_prep_date_gs']}.csv.gz")
        hgnc_ids_for_prediction = list(feat_test['hgnc_id_int'])
    
    # machine learning training set
    elif gene_set == 'train':
        feat_train = pd.read_csv(f"gene_score/training/train_test_data/feat_train_{config_vars['data_prep_date_gs']}.csv.gz")
        hgnc_ids_for_prediction = list(feat_train['hgnc_id_int'])
    
    # machine learning training and test set
    elif gene_set == 'train_test':
        feat_test = pd.read_csv(f"gene_score/training/train_test_data/feat_test_{config_vars['data_prep_date_gs']}.csv.gz")
        feat_train = pd.read_csv(f"gene_score/training/train_test_data/feat_train_{config_vars['data_prep_date_gs']}.csv.gz")
        hgnc_ids_for_prediction = list(feat_test['hgnc_id_int']) + list(feat_train['hgnc_id_int'])
    
    # all genes 
    elif gene_set == 'all':
        hgnc_ids_for_prediction = list(all_genes_scaled['hgnc_id_int'])
    
    # genes without known labels (all but training and test genes)
    elif gene_set == 'unknown':
        feat_test = pd.read_csv(f"gene_score/training/train_test_data/feat_test_{config_vars['data_prep_date_gs']}.csv.gz")
        feat_train = pd.read_csv(f"gene_score/training/train_test_data/feat_train_{config_vars['data_prep_date_gs']}.csv.gz")
        hgnc_ids_for_prediction = list(set(all_genes_scaled['hgnc_id_int']) - set(feat_test['hgnc_id_int']) - set(feat_train['hgnc_id_int']))
        
    # error in case of undefined/invalid gene set     
    else:
        raise ValueError("'gene_set' must be 'train', 'test', 'train_test', 'all' or 'unknown'.")
    
    # filter genes for prediction
    genes_for_predictions = all_genes_scaled.query("hgnc_id_int in @hgnc_ids_for_prediction")

    return genes_for_predictions



# In[18]:


def get_symbol_from_hgnc_id_int(hgnc_id_int_list):
    """
    Function that annotates a gene list of HGNC IDs with their symbol based 
    on the HGNC table of kidney-genetics on github.
    """
    # HGNC table gitHub raw URL 
    url = f"https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/B_AnnotationHGNC/results/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version_vs']}.csv.gz"

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


# In[19]:


def get_gene_set_from_hgnc_id_int(hgnc_id_int_list):
    """
    Function that annotates a gene list of HGNC IDs with their gene set (train, test, None).
    """
    test_ids = pd.read_csv(f"gene_score/training/train_test_data/feat_test_{config_vars['data_prep_date_gs']}.csv.gz")['hgnc_id_int'].tolist()
    train_ids = pd.read_csv(f"gene_score/training/train_test_data/feat_train_{config_vars['data_prep_date_gs']}.csv.gz")['hgnc_id_int'].tolist()
    
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


# In[20]:


def get_evidence_count_from_hgnc_id_int(hgnc_id_int_list):
    """
    Function that annotates a gene list of HGNC IDs with their evidence count.
    """
    # get positive and dispensible/negative genes
    pos_genes = pd.read_csv(f"gene_score/labels/results/positive_genes_{config_vars['creation_date_gs']}.csv.gz")
    neg_genes = pd.read_csv(f"gene_score/labels/results/dispensible_genes_{config_vars['creation_date_gs']}.csv.gz")
    
    # set evidence count of dispensible genes to -1
    neg_genes['evidence_count'] = -1

    # rowbind both
    annotated_genes = pd.concat([pos_genes[['hgnc_id_int', 'evidence_count']], neg_genes], ignore_index=True)
    
    # create df from gene list
    genes_df = pd.DataFrame({'hgnc_id_int' : hgnc_id_int_list})
    
    # annotate genes df with evidence count
    genes_df = genes_df.merge(annotated_genes[['hgnc_id_int', 'evidence_count']], how='left', left_on='hgnc_id_int', right_on='hgnc_id_int')
    
    return genes_df['evidence_count'].tolist()


# In[21]:


def make_predictions(ID, gene_set, save):
    # get configuration dictionary and results dictionary of the chosen experiment
    config_dic, results_dic = get_config_results_dics(ID=ID, score='gs') 
    
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
        NGS.to_csv(f"gene_score/predictions/results/NGS_predictions_ID{ID}_{gene_set}_{current_date}.csv.gz", index=False, compression='gzip')
        
    return NGS


# In[9]:


# make predictions for final chosen model for training set
NGS = make_predictions(ID=97, gene_set='train', save=False) 
NGS['ec2345'] = np.where(NGS['evidence_count'] == -1, 0, 1)

roc_auc_score(NGS['ec2345'], NGS['NGS'])


# In[20]:


# make predictions for final chosen model for test set
NGS = make_predictions(ID=97, gene_set='test', save=True) 
NGS['ec2345'] = np.where(NGS['evidence_count'] == -1, 0, 1)

# calculate AUC
roc_auc_score(NGS['ec2345'], NGS['NGS'])


# In[19]:


## Plot AUCs
from sklearn.metrics import roc_curve, auc


fpr1, tpr1, _ = roc_curve(NGS['ec2345'], NGS['NGS'])
roc_auc1 = auc(fpr1, tpr1)

# plotting
plt.figure()
plt.plot(fpr1, tpr1, color='darkorange', lw=2, label=f'{gene_set}(AUC = %0.2f)' % roc_auc1)

plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc='lower right')
plt.show()


# In[38]:


# Boxplot

gene_set = 'all' # TODO: add as variable in function!!!
ID = 97 # TODO: add as variable in function!!!

# NGS = make_predictions(ID=ID, gene_set='all', save=True) # TODO: add as variable in function!!!
NGS = make_predictions(ID=ID, gene_set=gene_set, save=False) # TODO: add as variable in function!!!



## create NCS boxplots based on evidence counts
evidence_counts = [-1, 0, 1, 2, 3, 4, 5]

plt.figure(figsize=(7, 6))
for count in evidence_counts:
    subset_df = NGS[NGS['evidence_count'] == count]
    plt.boxplot(subset_df['NGS'], positions=[count], widths=0.6, showfliers=True, vert=True)

    # Add count above the boxplot
    num_genes = len(subset_df)
    if num_genes > 0:
        plt.text(count, -0.05, f'{num_genes}', ha='center', va='bottom')
        
# add the boxplot for evidence_count = nan
subset_nan = NGS[NGS['evidence_count'].isna()]
if not subset_nan.empty:
    boxplot_nan = plt.boxplot(subset_nan['NGS'], positions=[6], widths=0.6, showfliers=True, vert=True, boxprops=dict(color='red'))
    num_genes_nan = len(subset_nan)
    plt.text(6, -0.05, f'{num_genes_nan}', ha='center', va='bottom')
        
# plot the gene number label
plt.text(7.05, -0.05, 'no. genes', ha='center', va='bottom')

# add a horizontal line
plt.axhline(y=-0.01, color='black', linestyle='--', xmin=0, xmax=1)

        
# Set labels and title
plt.xlabel('Evidence Count$^*$', fontsize=14)
plt.ylabel('Nephro Gene Score', fontsize=14)
# plt.title('Boxplots of Nephro Gene Score by Evidence Count', fontsize=15)

# Customize x-axis tick labels
custom_labels = {6: 'novel \n genes', -1: 'dispensable \n genes'}  # Set 'apple' label for x-position 7
tick_labels = [custom_labels.get(pos, pos) for pos in [-1, 0, 1, 2, 3, 4, 5, 6]]
plt.xticks([-1, 0, 1, 2, 3, 4, 5, 6], tick_labels, fontsize=12)
plt.yticks(fontsize=12)

# save plot
current_date = datetime.now().strftime("%Y-%m-%d")
plt.savefig(f"gene_score/predictions/results/boxplots_NGS_by_EC_ID{ID}_{gene_set}_{current_date}.png", dpi = 450, format='png')
plt.savefig(f"gene_score/predictions/results/boxplots_NGS_by_EC_ID{ID}_{gene_set}_{current_date}.pdf", dpi = 450, format='pdf')


# Show the plot
plt.show()


# In[77]:


# write NGS to json file per gene
# ID = 97
# convert to JSON with column names as keys and values in a list
# NGS = make_predictions(ID=ID, gene_set='all', save=True) # TODO: add as variable in function!!!

NGS_filepath = "predictions/results/NGS_predictions_ID97_all_2024-03-22.csv.gz"

NGS_for_json = pd.read_csv(NGS_filepath)

# replace NaN values of 'gene_set' with "none"
NGS_for_json['gene_set'] = NGS_for_json['gene_set'].fillna('none')

# replace NaN values of 'evidence_count' with None
NGS_for_json['evidence_count'] = NGS_for_json['evidence_count'].replace({np.nan: None})

# rename df to camelCase
NGS_for_json = NGS_for_json.rename(columns={'hgnc_id_int': 'hgncIdInt',
                                            'NGS':'ngs',
                                            'symbol': 'symbol',
                                            'gene_set': 'geneSet',
                                            'evidence_count': 'evidenceCount'})

# add metadata
meta_date = NGS_filepath.split("_")[-1].split(".")[0]
meta_ID = int([i for i in NGS_filepath.split("_") if i.startswith("ID")][0].split("ID")[-1]) # TODO: change??
config_dic, results_dic = get_config_results_dics(ID=int(meta_ID))
meta_clf = str(type(config_dic['clf']).__name__)
meta_version = "0.1.0" # TODO: specifiy manually??

meta_dic = {'date': meta_date, 'classifier': meta_clf, 'ID': meta_ID, 'version': meta_version}

for row in np.arange(NGS_for_json.shape[0]):
    json_data = NGS_for_json.iloc[row].to_dict()
    json_data['meta'] = meta_dic
    hgnc_id_int = json_data['hgncIdInt']
    symbol = json_data['symbol']

    # save to a file
    with open(f'predictions/results/json/hgnc_new/{hgnc_id_int}.json', 'w') as json_file:
        json.dump(json_data, json_file)
        
    # save to a file
    with open(f'predictions/results/json/symbols_new/{symbol}.json', 'w') as json_file:
        json.dump(json_data, json_file)


# In[9]:


# create an index for hgnc ids

# Directory containing the JSON files
folder_path = 'gene_score/predictions/results/json/hgnc/'

# List to store numbers extracted from file names
numbers = []

# Iterate over files in the directory
for filename in os.listdir(folder_path):
    # Check if file is JSON
    if filename.endswith('.json'):
        # Extract number from filename
#         number = int(filename.split('.')[0]) # Extract number as numeric
        number = filename.split('.')[0]  # Extract number as string

        numbers.append(number)

# Sort the numbers
# numbers.sort()
numbers.sort(key=int) # if numbers are extracted as string

# Create the hgnc_index.json file
output_file = 'gene_score/predictions/results/json/hgnc_index.json'
with open(output_file, 'w') as f:
    json.dump(numbers, f, indent=4)

print(f'Successfully created {output_file} with sorted numbers: {numbers}')


# In[10]:


# create summary json file with info of all genes

# Directory containing the JSON files
folder_path = 'gene_score/predictions/results/json/symbols/'

# List to store all gene data
gene_data_list = []

# Iterate over files in the directory
for filename in os.listdir(folder_path):
    if filename.endswith('.json'):
        file_path = os.path.join(folder_path, filename)
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                gene_data_list.append(data)
        except json.JSONDecodeError as e:
            print(f"Error reading {filename}: {e}")

# Create the summary JSON file
output_file = 'gene_score/predictions/results/json/gene_info_summary.json'
with open(output_file, 'w') as f:
    json.dump(gene_data_list, f, indent=4)

print(f'Successfully created {output_file} with {len(gene_data_list)} entries.')


# In[37]:


# Plot N-GS distribution of novel genes
a = pd.read_csv("gene_score/predictions/results/NGS_predictions_ID97_all_2024-03-22.csv.gz")
a = a.query("evidence_count.isna()")

plt.hist(a['NGS'], bins=100, density=True)
plt.xlabel('Nephro Gene Score', fontsize=14)
plt.ylabel('No. genes', fontsize=14)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.title('Distribution of Nephro Gene Score \n for Novel Genes', fontsize=14)
plt.savefig(f"gene_score/predictions/results/distribution_NGS_novel_genes_ID{ID}_{gene_set}_{current_date}.png", dpi = 450, format='png')
plt.savefig(f"gene_score/predictions/results/distribution_NGS_novel_genes_ID{ID}_{gene_set}_{current_date}.pdf", dpi = 450, format='pdf')

