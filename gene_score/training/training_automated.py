#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import config
from config_ML import *

# import basic modules
import numpy as np
import pandas as pd

# set options
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

# import preprocessing functions
from helper_functions_ML import *

# import modules for training
# from sklearn.model_selection import GridSearchCV
# from sklearn.metrics import roc_auc_score

# import classifiers
from sklearn.ensemble import AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.naive_bayes import GaussianNB


# In[2]:


print(get_latest_ML_config_training_file(config_dir))
pd.read_csv(get_latest_ML_config_training_file(config_dir))


# In[ ]:


# set configurations for ML experiment


# # 
# estimator = None
# clf = GaussianNB()
# model = "NaiveBayes"
# param_grid = {
#     'var_smoothing': [1e-12, 1e-9, 1e-5], # first parameter displayed on the y-axis in the 2D heatmap
#     'priors': [None], # second paramteter displayed on the x-axis in the 2D heatmap
# }


estimator = None
clf = SVC()
model = "SVM"
# Define a parameter grid for GridSearch
param_grid = {
    'C': np.logspace(-5, 5, 11),
    'gamma': np.logspace(-5, 5, 11),
    'kernel': ['poly'],
    'degree': [3]
}


# estimator = None
# clf = KNeighborsClassifier(n_neighbors=3)
# model = "KNN"
# param_grid = {
#     'n_neighbors': np.arange(10, 100),
#     'weights': ['uniform', 'distance'],
# }


# estimator = None
# clf = LinearSVC(dual=False)
# model = "SVM"
# param_grid = {
#     'C': np.logspace(-4, 3, 20), # first parameter displayed on the y-axis in the 2D heatmap
#     'penalty': ['l1', 'l2'], # second paramteter displayed on the x-axis in the 2D heatmap
#     'loss': ['squared_hinge']
# }



# estimator = DecisionTreeClassifier(max_depth=3)
# clf = AdaBoostClassifier(estimator=estimator)
# param_grid = {
#     'n_estimators': [5, 10, 50],
#     'learning_rate': np.linspace(0.01, 0.5, 10)
# }

# estimator = None
# clf = SVC()
# model = "SVM"
# param_grid = {
#     'C': np.logspace(2, 4, 13),
#     'kernel': ['rbf'],
#     'gamma': np.logspace(-7, -4, 15)
# }
cv = 5
scoring = 'roc_auc'
feature_groups_selected = ['gnomad', 'cellxgene', 'descartes', 'gtex', 'mgi', 'paralogues', 'phasCons', 'CpG_o2e']
drop_features = []
omit_scaling_features = ['paralogues', 'mgi']
scaling = 'robust'
pca_components = False
additional_info = '' 

# create config files
ID = create_ML_config(config_dir = 'config_files',
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


# # plot grid with performance score
# plot_2D_heatmap(ID=config_dic['ID'],
#                 cv_results=cv_results,
#                 param1=get_heatmap_params(param_grid=config_dic['param_grid'])[0],
#                 param2=get_heatmap_params(param_grid=config_dic['param_grid'])[1],
#                 figsize=(len(config_dic['param_grid'][get_heatmap_params(param_grid=config_dic['param_grid'])[0]]) + 5, 
#                          len(config_dic['param_grid'][get_heatmap_params(param_grid=config_dic['param_grid'])[1]]) + 3),
#                 save=True,
#                 show=True
#                )

# # get feature importance by permutation importance
# feature_imp = get_permutation_importance(ID=config_dic['ID'],
#                                          X_train = X_train, 
#                                          y_train = y_train, 
#                                          features = config_dic['features'],
#                            classifier = best_classifier,
#                            scoring=config_dic['scoring'],
#                            plot=True)



# In[10]:


cv_results


# In[14]:


plot_2D_heatmap(ID=config_dic['ID'],
                cv_results=cv_results,
                param1=get_heatmap_params(param_grid=config_dic['param_grid'])[0],
                param2=get_heatmap_params(param_grid=config_dic['param_grid'])[1],
                figsize=( 20, 
                         50),
                save=True,
                show=True
               )


# In[18]:


# get feature importance by permutation importance
feature_imp = get_permutation_importance(ID=config_dic['ID'],
                                         X_train = X_train, 
                                         y_train = y_train, 
                                         features = config_dic['features'],
                           classifier = best_classifier,
                           scoring=config_dic['scoring'],
                           plot=True)


# In[ ]:




