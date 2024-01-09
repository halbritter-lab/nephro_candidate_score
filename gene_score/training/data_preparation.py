# import basic modules
import numpy as np
import pandas as pd
import yaml
import os

# import preprocessing modules
from sklearn.model_selection import train_test_split

# get config file
CONFIG_FILE = "config_NCS.yml"   #TODO: CONFIG_FILE = os.getenv('CONFIG_FILE')

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"
script_path = "/gene_score/training/"

# read configs
with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}{script_path}")

# load features and labels
raw_feat = pd.read_csv(f"../features/results/gene_features_{config_vars['creation_date']}.csv.gz")
pos_genes = pd.read_csv(f"../labels/results/positive_genes_{config_vars['creation_date']}.csv.gz")
neg_genes = pd.read_csv(f"../labels/results/dispensible_genes_{config_vars['creation_date']}.csv.gz")

# load table with filling methods for missing values
filling_methods = pd.read_csv(f"raw/fill_missing_methods_{config_vars['creation_date']}.csv.gz", delimiter=";")

# create new column indicating if a evidence_count > 1
pos_genes['ec2345'] = np.where(pos_genes['evidence_count'] > 1, 1, 0)

# filter only positive genes with evidence_count > 1
pos_genes = pos_genes.query("ec2345 == 1")

# set negative genes to ec2345 = 0 
neg_genes['ec2345'] = 0

# combine positive and negative genes
labels = pd.concat([neg_genes[['hgnc_id_int', 'ec2345']], pos_genes[['hgnc_id_int', 'ec2345']]])

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
labels_train.to_csv(f"train_test_data/labels_train_{config_vars['data_prep_date']}.csv.gz", index=False, compression='gzip')
labels_test.to_csv(f"train_test_data/labels_test_{config_vars['data_prep_date']}.csv.gz", index=False, compression='gzip')
feat_train.to_csv(f"train_test_data/feat_train_{config_vars['data_prep_date']}.csv.gz", index=False, compression='gzip')
feat_test.to_csv(f"train_test_data/feat_test_{config_vars['data_prep_date']}.csv.gz", index=False, compression='gzip')

# create table with median values of training set for filling NA
median_values_train = feat_train.drop(columns = ['hgnc_id_int']).median().reset_index()
median_values_train.columns = ['feature', 'median']

# save median df
median_values_train.to_csv(f"train_test_data/median_values_train_{config_vars['data_prep_date']}.csv.gz", index=False, compression='gzip')
