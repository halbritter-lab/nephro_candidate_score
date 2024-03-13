# import basic modules
import sys
import numpy as np
import pandas as pd
import yaml

# import preprocessing functions
from helper_functions_ML import *

# import classifiers
from sklearn.ensemble import AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.linear_model import RidgeClassifier
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier


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


def main():
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python get_SHAP_values.py --ID <value>")
        sys.exit(1)

    # Extract the value of myID
    if sys.argv[1] == "--ID":
        ID = int(sys.argv[2])
    else:
        print("Invalid argument. Use --ID.")
        sys.exit(1)
        
    config_dic, results_dic = get_config_results_dics(ID=ID)  
    
    
    get_SHAP_values(ID=ID, 
                    X=config_dic['X_train'], 
                    y=config_dic['y_train'], 
                    features=config_dic['features'], 
                    clf=results_dic['best_classifier'], 
                    create_plot_pdf=True)
    
    
    
    
#     # run permutation importance five times and rank features 
#     df_combined = pd.DataFrame()

#     for i in np.arange(5): 
#         print(f"Getting permutation importance with random_state = {i}")
#         perm_imp = get_permutation_importance(ID=ID,
#                                               X_train=config_dic['X_train'], 
#                                               y_train=config_dic['y_train'], 
#                                               features=config_dic['features'], 
#                                               classifier=results_dic['best_classifier'],
#                                               scoring=config_dic['scoring'],
#                                               plot=False,
#                                               random_state=i
#                                              )
    
    
#         perm_imp['rank_value'] = 1 / perm_imp['perm_imp'].rank(ascending=False)

#         df_combined = pd.concat([df_combined, perm_imp], ignore_index=True)

#     sum_rank_values = df_combined.groupby('feature')['rank_value'].sum().reset_index().sort_values(by='rank_value', ascending=False)
    
    


if __name__ == "__main__":
    main()




