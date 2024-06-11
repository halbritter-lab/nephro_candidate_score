# import basic modules
import numpy as np
import pandas as pd
import sys
import yaml

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

# append path where common functions are located
sys.path.append(f"{config_vars['ML_projectsdir']}{project_name}")

# import common functions
from common_functions.training_helper_functions import *

# get config file
TODO: CONFIG_FILE = os.getenv('CONFIG_FILE')

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"

# read configs
with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}")


# set score type
score = 'gs' # TODO: customize
score_string, id_string, label_string = get_score_specific_args(score)


# define main function 
def main():
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python get_permutation_importance.py --ID <value>")
        sys.exit(1)

    # Extract the value of myID
    if sys.argv[1] == "--ID":
        ID = int(sys.argv[2])
    else:
        print("Invalid argument. Use --ID.")
        sys.exit(1)
        
    config_dic, results_dic = get_config_results_dics(ID=ID)  
    
    # run permutation importance five times and rank features 
    df_combined = pd.DataFrame()

    for i in np.arange(5): 
        print(f"Getting permutation importance with random_state = {i}")
        perm_imp = get_permutation_importance(
            ID=ID,
            X_train=config_dic['X_train'], 
            y_train=config_dic['y_train'], 
            features=config_dic['features'], 
            classifier=results_dic['best_classifier'],
            scoring=config_dic['scoring'],
            plot=False,
            random_state=i,
            score=score
        )
    
    
        perm_imp['rank_value'] = 1 / perm_imp['perm_imp'].rank(ascending=False)

        df_combined = pd.concat([df_combined, perm_imp], ignore_index=True)

    sum_rank_values = df_combined.groupby('feature')['rank_value'].sum().reset_index().sort_values(by='rank_value', ascending=False)
    
    
    # save results
    date_time = datetime.today().strftime('%Y-%m-%d--%H-%M-%S')
    sum_rank_values.to_csv(f"{score_string}/training/feature_importance/permutation_importance/perm_imp_sum_rank_ID{ID}_{date_time}.csv", index=False)


if __name__ == "__main__":
    main()




