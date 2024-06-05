# import basic modules
import numpy as np
import os
import pandas as pd
import sys
import yaml

# import sklearn classifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"
# script_path = "/variant_score/"

# read configs
CONFIG_FILE = os.getenv('CONFIG_FILE')

with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
# os.chdir(f"{config_vars['ML_projectsdir']}{project_name}{script_path}")
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}")

# append path where common functions are located
sys.path.append(f"{config_vars['ML_projectsdir']}{project_name}")


# import common functions
# from helper_functions_ML_vs import *
from common_functions.training_helper_functions import *

score = 'vs'
score_string, id_string, label_string = get_score_specific_args(score)

# set smote
smote = False

# # # set classifer and param grid for training - Decision Tree 
# estimator = None
# clf = DecisionTreeClassifier(random_state=1) 
# model = "DecisionTree"
# param_grid = {
#     'criterion' : ['gini', 'entropy', 'log_loss'], 
#     'max_depth': [3, 4, 5, 6, 7, 8, 9, 10], 
#     'min_samples_split': [30, 60, 90, 120, 150],
#     'min_samples_leaf' : [4,5,6,7,8,9,10]
# #     'max_leaf_nodes' :[]
# }


estimator = None
clf = LogisticRegression(random_state=1, max_iter=1000) 
model = "logReg"
param_grid = {
    'penalty' : [None, 'l2', 'l1', 'elasticnet'], 
#     'C': [100, 1000, 10000, 100000, 1000000, 10000000, 100000000]
    'C':[100]
    }



IMPACT_prop=True
additional_info = 'IMPACT_prop' 
additional_info = 'delete' 




# set names of selected feature groups
feature_groups_selected = ['Consequence',
#                            'NMD',
                           'gnomAD_AF',
#                            'SpliceAI',
                           'CADD',
                           'IMPACT',
#                            'phastCons'
                          ]





# set configurations
cv = 5
scoring = 'roc_auc'
drop_features = ['Csq_3_prime_UTR_variant', 
                 'Csq_5_prime_UTR_variant', 
                 'Csq_NMD_transcript_variant', 
                 'Csq_coding_sequence_variant', 
                 'Csq_downstream_gene_variant', 
#                  'Csq_frameshift_variant', 
                 'Csq_incomplete_terminal_codon_variant', 
                 'Csq_inframe_deletion', 
                 'Csq_inframe_insertion', 
                 'Csq_intron_variant', 
#                  'Csq_missense_variant', 
                 'Csq_non_coding_transcript_exon_variant', 
                 'Csq_non_coding_transcript_variant', 
                 'Csq_protein_altering_variant', 
                 'Csq_splice_acceptor_variant', 
                 'Csq_splice_donor_5th_base_variant', 
                 'Csq_splice_donor_region_variant', 
                 'Csq_splice_donor_variant', 
                 'Csq_splice_polypyrimidine_tract_variant', 
                 'Csq_splice_region_variant', 
                 'Csq_start_lost', 
                 'Csq_start_retained_variant', 
#                  'Csq_stop_gained', 
                 'Csq_stop_lost', 
                 'Csq_stop_retained_variant', 
#                  'Csq_synonymous_variant', 
                 'Csq_upstream_gene_variant'
                ] # ['gnomADg_AF']
# omit_scaling_features = ['Consequence', 'NMD', 'IMPACT']
# omit_scaling_features = ['IMPACT']
omit_scaling_features = []


scaling = None
scaling = 'robust'
scaling = 'standard'


pca_components = False



# create config files
ID = create_ML_config(config_dir = f"{score_string}/training/config_files/",
                     estimator = estimator,
                     clf = clf, 
                     param_grid = param_grid, 
                     cv = cv, 
                     scoring = scoring,
                     feature_groups_selected = feature_groups_selected,
                     drop_features = drop_features,
                     omit_scaling_features = omit_scaling_features,
                     model = model, # defines which technique should be used for filling NaN
                     scaling = scaling,
                     pca_components = pca_components,
                     additional_info = additional_info,
                     score = score
                    )

print(f"Training_ID: ID{ID}")

# get config_dic
with open(f"{score_string}/training/config_files/config_dic_ID{ID}.pkl", 'rb') as file:
    config_dic = pickle.load(file)


# get the training data
X_train, y_train, features = get_training_data(model=config_dic['model'],
                                               feature_groups_selected=config_dic['feature_groups_selected'],
                                               drop_features=config_dic['drop_features'],
                                               omit_scaling_features=config_dic['omit_scaling_features'],
                                               scaling=config_dic['scaling'],
                                               pca_components=config_dic['pca_components'],
                                               IMPACT_prop=IMPACT_prop,
                                               score=score
                                              )

# fill config_dic with training data and features
config_dic['features'] = features
config_dic['X_train'] = X_train
config_dic['y_train'] = y_train

# dump config_dic
with open(f'{score_string}/training/config_files/config_dic_ID{ID}.pkl', 'wb') as file:
    pickle.dump(config_dic, file)

    
print("ready for training")    
# train model
cv_results, best_params, best_classifier = train_with_grid_search(ID=config_dic['ID'],
                                                                  X_train=X_train, 
                                                                  y_train=y_train, 
                                                                  estimator=config_dic['clf'], 
                                                                  param_grid=config_dic['param_grid'], 
                                                                  cv=config_dic['cv'], 
                                                                  scoring=config_dic['scoring'], 
                                                                  verbose=2,
                                                                  smote=smote,
                                                                  score=score
                                                                 )
