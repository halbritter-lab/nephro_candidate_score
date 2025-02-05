# import config

# import basic modules
import sys
import numpy as np
import os
import pandas as pd
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
from xgboost import XGBClassifier

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"

# read configs
CONFIG_FILE = os.getenv('CONFIG_FILE')

with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}")

# append path where common functions are located
sys.path.append(f"{config_vars['ML_projectsdir']}{project_name}")

# import common functions
from common_functions.training_helper_functions import *

score = 'gs'
score_string, id_string, label_string = get_score_specific_args(score)


# set options
drop_features = False
pca_reduced_features = True
smote = False


# set classifer and param grid for training

# ## XGBoost
# estimator = None
# clf = XGBClassifier(random_state=1, booster='gbtree')
# model = 'DecisionTree'
# param_grid = {
#     'n_estimators': [50, 70, 90], #np.arange(70, 180, 30),
#     'max_depth' : [3,4,5,6], #np.arange(1, 5, 1),
#     'learning_rate': np.linspace(1e-1, 2.5e-1, 3), #np.logspace(-2, -1, 5),
#     'reg_alpha' : np.linspace(3.2e-1, 5, 4), #np.logspace(-2, 1, 5),
#     'reg_lambda' : np.logspace(-4, -2, 4),
#     'subsample': [0.8, 0.85, 0.9],
#     'gamma': np.logspace(-2, 1, 4)
# }

   
# ## AdaBoost classifier with estimator=DecisionTree
# estimator = DecisionTreeClassifier(max_depth=2) #TODO: adapt max_depth
# clf = AdaBoostClassifier(estimator=estimator, random_state=1)
# model = "DecisionTree"
# param_grid = {
#     'n_estimators': np.arange(5, 50, 5),
#     'learning_rate': np.logspace(-2, 0, 10)
# }


# ## Decision Tree
# estimator = None
# clf = DecisionTreeClassifier(random_state=1) 
# model = "DecisionTree"
# param_grid = {
#     'max_depth': [3], #[3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 
#     'min_samples_split': [40], #[40, 50, 60, 70, 80, 90, 100, 110],
#     'min_samples_leaf' : [3] #[3,4,5,6,7]
# #     'max_leaf_nodes' :[]
# }


# ## Random Forest
# estimator = None
# clf = RandomForestClassifier()
# model = "DecisionTree"
# param_grid = {
#     'n_estimators': np.arange(50, 250, 25),
#     'max_depth': np.arange(10, 50, 5),
#     'min_samples_leaf': np.arange(2, 10, 1),
#     'min_samples_split': np.arange(2, 50, 2)
# }


# ## Ridge Classifier
# estimator = None
# clf = RidgeClassifier()
# model = "logReg"
# param_grid = {
#     'alpha': np.logspace(3, 5, 10), 
#     'tol': [1e-6,  1e-5]
# }


# ## Gaussian Process Classifier
# estimator = None
# clf = GaussianProcessClassifier(random_state=1)
# model = "GPC"
# param_grid = {
#     'optimizer': ['fmin_l_bfgs_b', None], #[ 0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10, 100],
#     'n_restarts_optimizer': [0, 5, 10, 15]
# }


# ## Quadratic Discriminant Analysis
# estimator = None
# clf = QuadraticDiscriminantAnalysis()
# model = "QDA"
# param_grid = {
#     'reg_param': np.logspace(-2, 0, 10), #[ 0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10, 100],
#     'tol': [1e-2, 0, 1e+2 ]
# }


# ## MLP Classifier
# estimator = None
# clf = MLPClassifier(random_state=1)
# model = "MLP"
# param_grid = {
#     'alpha': np.logspace(-2, 3, 10),
#     'hidden_layer_sizes': [(70,), (60,), (50,), (40,), (50, 10), (50, 20)],
#     'solver': ['lbfgs'],
#     'activation': ['logistic', 'tanh', 'relu']
# }


# ## KNN
# estimator = None
# clf = KNeighborsClassifier()
# model = "KNN"
# param_grid = {
#     'n_neighbors': np.arange(10, 100, 2), 
#     'weights' : ['uniform', 'distance']
# }


# ## SVM - linear
# estimator = None
# clf = LinearSVC(dual=False)
# model = "SVM"
# param_grid = {
#     'C': np.logspace(-4, 3, 20),
#     'penalty' : ['l1', 'l2'],
#     'loss': ['squared_hinge']
# }


# ## SVM - rbf
# estimator = None
# clf = SVC()
# model = "SVM"
# param_grid = {
#     'C': np.logspace(2, 4, 13),
#     'kernel': ['rbf'],
#     'gamma': np.logspace(-7, -4, 15)
# }


# ## SVM - poly
# estimator = None
# clf = SVC(random_state=1)
# model = "SVM"
# param_grid = {
#     'C': np.logspace(-2, 8, 9),
#     'gamma': np.logspace(-5, -0, 9),
#     'kernel': ['poly'],
#     'degree': [3]
# }

estimator = None

model = 'SVM'
clf = SVC(random_state=1)

param_grid = {
    'C': np.logspace(-2, 8, 9),
    'gamma': np.logspace(-5, -0, 9),
    'kernel': ['poly'],
    'degree': [3]
}


cv = 5
scoring = 'roc_auc'
feature_groups_selected = ['gnomad', 'cellxgene', 'descartes', 'gtex', 'mgi', 'paralogues', 'phasCons', 'CpG_o2e']
drop_features = []
omit_scaling_features = ['paralogues', 'mgi']
scaling = 'robust'
pca_components = False
additional_info = '' 
additional_info = 'delete' # TODO: remove



###################################################################################################### 
##### OPTION: repeat specific experiment with only most important features (drop other features) #####

if drop_features:

    # specify experiment ID of which most important features should be taken for new experiment  
    previous_ID = 41 

    # specify how many importamt features should be used
    index = 25

    # load configuration dictionary of previous experiment
    config_dic, results_dic = get_config_results_dics(ID=previous_ID, score=score) 

    # select drop features
    important_feat, unimportant_feat = select_feat_from_perm_imp(ID=previous_ID, index=index, score=score)    
    drop_features = unimportant_feat

    # use configuarations as in the previous experiment
    model = config_dic['model']
    param_grid = config_dic['param_grid']
    cv = config_dic['cv']
    scoring = config_dic['scoring']
    feature_groups_selected = config_dic['feature_groups_selected']
    omit_scaling_features = config_dic['omit_scaling_features']
    scaling = config_dic['scaling']
    pca_components = False
    additional_info = ''
    additional_info = 'delete' # TODO: remove
######################################################################################################


###################################################################################################### 
##### OPTION: train with PCA reduced feature set #####
# Note: the PCA is not performed on all features, but a separate PCA is performed on feature groups, 
# with high correlation, e.g. gnomad features. The respective PCA reduced feature set was calculated 
# before and is loaded in get_feat_reduced_trainig_data(...)

if pca_reduced_features:

    # set percentage of variance explained by principal components until which PC should be kept
    perc_var_exp = 95

    # set names of selecte feature groups (needed for indication in config_ML_training_XXX.csv to indicate 
    # that PCA reduced feature groups where used)
    feature_groups_selected = [f'gnomad_PC_{perc_var_exp}', 
                               f'cellxgene_PC_{perc_var_exp}', 
                               f'descartes_PC_{perc_var_exp}', 
                               f'gtex_PC_{perc_var_exp}', 
                               'mgi', 
                               'paralogues', 
                               'phasCons', 
                               'CpG_o2e']

######################################################################################################


# create config files
ID = create_ML_config(
    config_dir = f"{score_string}/training/config_files/",
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
    score=score
)

print(f"Training_ID: ID{ID}")

# get config_dic
with open(f"{score_string}/training/config_files/config_dic_ID{ID}.pkl", 'rb') as file:
    config_dic = pickle.load(file)



if pca_reduced_features:
    # get the training data
    X_train, y_train, features = get_feat_reduced_trainig_data(perc_var_exp=perc_var_exp,    
                                                               model=model, 
                                                               omit_scaling_features=omit_scaling_features,
                                                               scaling=scaling,
                                                              score=score)


else:
    # get the training data
    X_train, y_train, features = get_training_data(
        model=config_dic['model'],
        feature_groups_selected=config_dic['feature_groups_selected'],
        drop_features=config_dic['drop_features'],
        omit_scaling_features=config_dic['omit_scaling_features'],
        scaling=config_dic['scaling'],
        pca_components=config_dic['pca_components'],
        IMPACT_prop = False, # only needed vor variant score training
        score=score
)

# fill config_dic with training data and features
config_dic['features'] = features
config_dic['X_train'] = X_train
config_dic['y_train'] = y_train

# dump config_dic
with open(f"{score_string}/training/config_files/config_dic_ID{ID}.pkl", 'wb') as file:
    pickle.dump(config_dic, file)

    
# print("ready for training")    

# train model
cv_results, best_params, best_classifier = train_with_grid_search(
    ID=config_dic['ID'],
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




