# import config
from config_ML import *

# import basic modules
import sys
import numpy as np
import pandas as pd

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


# define ID, index of important features, classifier
ID = 41
index = 25

# load configuration dictionary for specified ID
config_dic, results_dic = get_config_results_dics(ID=ID)    

# set classifer and param grid for training
# estimator = config_dic['estimator']
model = config_dic['model']
# clf = RandomForestClassifier(random_state=1)
# clf = MLPClassifier(random_state=1)

# estimator = DecisionTreeClassifier(max_depth=2) #TODO: adapt max_depth
# clf = AdaBoostClassifier(estimator=estimator, random_state=1)

estimator = None
# clf = DecisionTreeClassifier(random_state=1)
# clf = SVC(random_state=1)
# clf = LinearSVC(dual=False, random_state=1)
# clf = KNeighborsClassifier()
# clf = QuadraticDiscriminantAnalysis()
# clf = RidgeClassifier(random_state=1)
clf = SVC(random_state=1)




# use configuarations as in the previous training
param_grid = config_dic['param_grid']
cv = config_dic['cv']
scoring = config_dic['scoring']
feature_groups_selected = config_dic['feature_groups_selected']
omit_scaling_features = config_dic['omit_scaling_features']
scaling = config_dic['scaling']
pca_components = False
additional_info = '' 

# select drop features
important_feat, unimportant_feat = select_feat_from_perm_imp(ID=ID, index=index)    
drop_features = unimportant_feat





# ## XGBoost
# estimator = None
# clf = XGBClassifier(random_state=1, booster='gbtree')
# model = 'DecisionTree'
# param_grid = {
#     'n_estimators': np.arange(10, 180, 30),
#     'max_depth' : np.arange(4, 9, 1),
#     'learning_rate': np.logspace(-2, 0, 5),
#     'reg_alpha' : np.logspace(-2, 1, 5),
#     'reg_lambda' : np.logspace(-4, -2, 5),
#     'subsample': [0.8, 0.9, 1.0],
#     'gamma': np.logspace(-3, 0, 5)
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
# clf = DecisionTreeClassifier() # TODO: set random_state =1
# model = "DecisionTree"
# param_grid = {
#     'max_depth': [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 
#     'min_samples_split': [40, 50, 60, 70, 80, 90, 100, 110],
#     'min_samples_leaf' : [3,4,5,6,7]
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
# clf = RidgeClassifier() #TODO: set random_state
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
# clf = SVC() # TODO: set random_state
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


# cv = 5
# scoring = 'roc_auc'
# feature_groups_selected = ['gnomad', 'cellxgene', 'descartes', 'gtex', 'mgi', 'paralogues', 'phasCons', 'CpG_o2e']
# drop_features = []
# omit_scaling_features = ['paralogues', 'mgi']
# scaling = 'robust'
# pca_components = False
# additional_info = '' 

# create config files
ID = create_ML_config(config_dir = config_dir,
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
                     additional_info = ''
                    )

print(f"Training_ID: ID{ID}")

# get config_dic
with open(f"{config_dir}/config_dic_ID{ID}.pkl", 'rb') as file:
    config_dic = pickle.load(file)


# get the training data
X_train, y_train, features = get_training_data(model=config_dic['model'],
                                     feature_groups_selected=config_dic['feature_groups_selected'],
                                               drop_features=config_dic['drop_features'],
                                     omit_scaling_features=config_dic['omit_scaling_features'],
                                     scaling=config_dic['scaling'],
                                     pca_components=config_dic['pca_components']
                                    )

# fill config_dic with training data and features
config_dic['features'] = features
config_dic['X_train'] = X_train
config_dic['y_train'] = y_train

# dump config_dic
with open(f'{config_dir}/config_dic_ID{ID}.pkl', 'wb') as file:
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
                                                                  verbose=2
                                                                 )
