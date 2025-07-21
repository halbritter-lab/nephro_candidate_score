#### Get SHAP values ####
# creates SHAP values for model of given ID


# import basic modules
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import pandas as pd
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

# import preprocessing functions
from helper_functions_ML_vs import *

# SHAP values
import shap
import warnings
from sklearn.preprocessing import StandardScaler


warnings.filterwarnings("ignore", message="is_sparse is deprecated and will be removed in a future version.")


ID = 2 # TODO: set ID
config_dic, results_dic = get_config_results_dics(ID=ID) 
features = config_dic['features']

random.seed(42)

X = config_dic['X_train'] 

X_lab = pd.DataFrame(X, columns=features)

y = config_dic['y_train']
clf = results_dic['best_classifier']
clf.random_state = 1 

model = clf.fit(X_lab, y)

explainer = shap.Explainer(model.predict, X_lab)  

shap_values = explainer(X_lab)

# dump SHAP values in pickle
with open(f"variant_score/training/feature_importance/SHAP_values/values/SHAP_values_ID{ID}_{datetime.today().strftime('%Y-%m-%d')}.pkl", 'wb') as file:
    pickle.dump(shap_values, file)





