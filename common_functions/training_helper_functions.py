# import basic modules
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import time
import yaml

# import sklearn modues
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn.inspection import permutation_importance # for feature importance by permutation importance
from imblearn.pipeline import Pipeline # for SMOTE technique
from imblearn.over_sampling import SMOTE # # for SMOTE technique


import shap

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"
# score_path = "/gene_score/"

# read configs
# CONFIG_FILE = "./training/config_NCS.yml" # TODO: change
CONFIG_FILE = os.getenv('CONFIG_FILE')

with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}")


class Tee:
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)

    def flush(self):
        for f in self.files:
            f.flush()


def is_numeric(x):
    return (all(isinstance(i, (int, float, np.number)) for i in x))



def get_score_specific_args(score):
    # check if score is correctly specified
    if score not in ['gs', 'vs']:
        raise ValueError(f"Invalid argument: {score}. Allowed values are: 'gs' for gene score or 'vs' for variant score.")
    
    score_string = 'gene_score' if score == 'gs' else 'variant_score'
    id_string = 'hgnc_id_int' if score == 'gs' else 'var_ID'
    label_string = 'ec2345' if score == 'gs' else 'P_LP'

    return score_string, id_string, label_string


def init_config_ML_training_file(config_dir, date_time, score):
    """
    Initializes a config_ML_training file with the current datetime, ID=0, and NaN as values.
    """
   
    score_string, _, _ = get_score_specific_args(score)
    
    # create configuration dictionary
    config_dic = {'date_time': [date_time],
                  'ID': [0],
                  'estimator': [None],
                  'clf': [None],
                  'param_grid': [None],
                  'cv': [None],
                  'scoring': [None],
                  'feature_groups_selected': [None],
                  'drop_features': [None],
                  'omit_scaling_features': [None],
                  'model': [None],
                  'scaling':[None],
                  'pca_components': [False],
                  'additional_info': ''
                 }
    
    init_df = pd.DataFrame(config_dic)
        
    # create csv
    init_df.to_csv(f"{score_string}/training/config_files/config_ML_training_{date_time}.csv", index=False)
    
    # set file to read-only
    os.chmod(f"{score_string}/training/config_files/config_ML_training_{date_time}.csv", 0o444)
    time.sleep(2)
    
    print(f"'config_ML_training_{date_time}.csv' initialized.")


def get_latest_ML_config_training_file(config_dir, score):
    """
    Returns the latest config_ML_training_XXX.csv file. 
    If none available, one is initialized.
    """
    
    score_string, _, _ = get_score_specific_args(score)
        
    # list all files in the config directory
    files = os.listdir(config_dir)

    # filter the files that match the naming pattern
    config_files = [file for file in files if file.startswith("config_ML_training")]

    # isolate the datetime parts
    date_times = [x.split('_')[-1].replace('.csv', '') for x in config_files]

    # parse the string into a datetime object
    date_objects = [datetime.strptime(x, "%Y-%m-%d--%H-%M-%S") for x in date_times]

    # get the latest config_ML_trainig_XXX.csv file, if none available initialize one
    if len(date_objects) == 0:
        date_time = datetime.today().strftime('%Y-%m-%d--%H-%M-%S')
        init_config_ML_training_file(config_dir, date_time, score)
        return(f"{score_string}/training/config_files/config_ML_training_{date_time}.csv")
        
    else:
        latest_dt = max(date_objects).strftime('%Y-%m-%d--%H-%M-%S')
        return(f"{score_string}/training/config_files/config_ML_training_{latest_dt}.csv")





def create_ML_config(config_dir,
                     estimator,
                     clf, 
                     param_grid, 
                     cv, 
                     scoring,
                     feature_groups_selected,
                     drop_features,
                     omit_scaling_features,
                     model, # defines which technique should be used for filling NaN
                     scaling,
                     pca_components,
                     additional_info,
                     score
                    ):
    """
    - dumps the given configurations in pickle with a new ID
    - updates the config_ML_XXX.csv
    """
    
    score_string, _, _ = get_score_specific_args(score)

    # latest config ML training file
    lt_cf = get_latest_ML_config_training_file(config_dir=config_dir, score=score)   
    
    # create new experiment ID (highest ID so far + 1)
    highest_ID = pd.read_csv(lt_cf)['ID'].max()
    new_ID = highest_ID + 1
    
    # get current date_time
    date_time = datetime.today().strftime('%Y-%m-%d--%H-%M-%S')
    
    # create a configuration dictionary for the experiment
    config_dic = {'date_time': date_time,
                  'ID': new_ID,
                  'estimator': estimator,
                  'clf': clf,
                  'param_grid': param_grid,
                  'params': clf.get_params(),
                  'cv': cv,
                  'scoring': scoring,
                  'feature_groups_selected': feature_groups_selected,
                  'drop_features': drop_features,
                  'omit_scaling_features': omit_scaling_features,
                  'model': model,
                  'scaling': scaling,
                  'pca_components': pca_components,
                  'additional_info': additional_info
                 }
    
    # save the configuration dictionary with pickle
    with open(f'{score_string}/training/config_files/config_dic_ID{new_ID}.pkl', 'wb') as file:
        pickle.dump(config_dic, file)
    print(f"{score_string}/training/config_dic_ID{new_ID}.pkl dumped.")
    
    # format entry - param_grid
    param_grid_string = '|'.join(
        f"{key}=[{value[0]}, {value[-1]}, {len(value)}]" if is_numeric(value)
        else f"{key}={str(value)}" 
        for key, value in param_grid.items()
    )
    
    # create a new config_ML_trainig_XXX.csv file with the current settings
    new_row = pd.DataFrame({'date_time': [date_time],
                            'ID': [new_ID], 
                            'estimator': [estimator],
                            'clf': [str(type(clf).__name__)],
                            'param_grid': [param_grid_string],
                            'cv': [cv], 
                            'scoring': [scoring],
                            'feature_groups_selected': ['|'.join(feature_groups_selected)],
                            'drop_features': ['|'.join(drop_features)],
                            'omit_scaling_features': ['|'.join(omit_scaling_features)],
                            'model': [model],
                            'scaling': [scaling],
                            'pca_components': [pca_components],
                            'additional_info': [additional_info]
                           })
    
    old_config = pd.read_csv(lt_cf)
    
    new_config = pd.concat([old_config, new_row], ignore_index=True)
   
    new_config.to_csv(f"{score_string}/training/config_files/config_ML_training_{date_time}.csv", index=False)
    
    # set file to read-only
    os.chmod(f"{score_string}/training/config_files/config_ML_training_{date_time}.csv", 0o444)
    
    print(f"New config file: 'config_ML_training_{date_time}.csv' created.")
    
    return(new_ID)


# function to select only specific feature names based on feature group
def get_features_from_groups(groups, feature_df, score):
    """
    e.g. groups = ['gnomad', 'descartes']
    
    """
    score_string, _, _ = get_score_specific_args(score)
    
    
    res = []
    
    
    if score_string == 'gene_score':
        feat_dic = {"gnomad": [x for x in feature_df.columns if x.startswith("gnomad")], 
                    "cellxgene": [x for x in feature_df.columns if x.startswith("CL_")],
                    "descartes": [x for x in feature_df.columns if x.endswith(("perc_expr", "_nTPM", "fetal_kidney_tau"))],
                    "gtex": [x for x in feature_df.columns if x.startswith("gtex")],
                    "mgi": [x for x in feature_df.columns if x.startswith("max_mgi_kid")],
                    "paralogues": [x for x in feature_df.columns if x.startswith("no_paralogues")],
                    "phasCons": [x for x in feature_df.columns if x.startswith("avg_phasCons")],
                    "CpG_o2e": [x for x in feature_df.columns if x.endswith("CpG_o2e_ratio")]
                   }
    else: # for variant score
        feat_dic = {"Consequence":[x for x in feature_df.columns if x.startswith("Csq_")],
                    "NMD":[x for x in feature_df.columns if x.startswith("NMD")],
                    "gnomAD_AF":[x for x in feature_df.columns if x.startswith("gnomAD")],
                    "SpliceAI":[x for x in feature_df.columns if x.startswith("SpliceAI")],
                    "CADD":[x for x in feature_df.columns if x.startswith("CADD")],
                    "IMPACT":[x for x in feature_df.columns if x.startswith("IMPACT")],
                    "phastCons":[x for x in feature_df.columns if x.startswith("phastCons")]
                   }

    for i in groups:
        if i not in feat_dic.keys():
            raise ValueError(f"Invalid group: '{i}'") 
        res = res + feat_dic[i]    
    
    return(res)   


# function to return the pre-defined method of filling missing values for each model and feature
def get_filling_method(model, feature, score):
    score_string, _, _ = get_score_specific_args(score)
    
    filling_methods = pd.read_csv(f"{score_string}/training/raw/fill_missing_methods_{config_vars[f'creation_date_{score}']}.csv.gz")
    method = filling_methods.loc[filling_methods['feature'] == feature, model].values[0]
    return(method)


# function to fill missing values with given method
def fill_missing_vals(df, model, score):
    score_string, _, _ = get_score_specific_args(score)
    
    df_filled = df.copy()
  
    # median values of training data
    median_values_train = pd.read_csv(f"{score_string}/training/train_test_data/median_values_train_{config_vars[f'data_prep_date_{score}']}.csv.gz")
    
    cols_with_missing = df.columns[df.isna().any()].tolist()
    for col in cols_with_missing: 
        # get filling method for specific feature and model
        method = get_filling_method(model, col, score)
        
        if pd.isna(method):
            raise ValueError(f"No filling method found for feature column '{col}' and model '{model}'")
        
        if method == 'zero':
            df_filled.fillna({col: 0}, inplace=True)
    
            
        elif method == 'median':
            med_val = median_values_train.loc[median_values_train['feature'] == col, 'median'].values[0]
#             df[col] = \
#             df[col].fillna(med_val)
            df_filled.fillna({col: med_val}, inplace=True)

    return(df_filled)


# function to numpy arrays of filled and scaled training data
def get_training_data(model,
                      feature_groups_selected,
                      drop_features, 
                      omit_scaling_features,
                      scaling, 
                      pca_components,
                      IMPACT_prop,
                      score
                     ):
    
    # get score and ID string
    score_string, id_string, label_string = get_score_specific_args(score)
    
    # load raw training data
    feat_train = pd.read_csv(f"{score_string}/training/train_test_data/feat_train_{config_vars[f'data_prep_date_{score}']}.csv.gz")
    labels_train = pd.read_csv(f"{score_string}/training/train_test_data/labels_train_{config_vars[f'data_prep_date_{score}']}.csv.gz")
    
    # get training data with balanced IMPACT groups
    if score == 'vs':
        if IMPACT_prop:
            feat_train = pd.read_csv(f"{score_string}/training/train_test_data/feat_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz")
            labels_train = pd.read_csv(f"{score_string}/training/train_test_data/labels_train_IMPACT_prop_{config_vars['data_prep_date_vs']}.csv.gz")
    
    
    # get selected features
    features_from_groups = get_features_from_groups(feature_groups_selected, feat_train, score)
    features_selected = [x for x in features_from_groups if not x in drop_features]
    feat_train = feat_train[[id_string] + features_selected]
    
    # fill features missing values
    feat_train_filled = fill_missing_vals(feat_train, model, score)
    
    # get features that should be scaled
    omit_scaling = get_features_from_groups(omit_scaling_features, feat_train, score) # features that should not be scaled
    scaling_features = [i for i in features_selected if i not in omit_scaling] # features that should be scaled
    
    if scaling == 'standard':
        # create StandardScaler
        stand_scal = StandardScaler()

        # scale features
        feat_train_scaled = feat_train_filled.copy()
        feat_train_scaled[scaling_features] = stand_scal.fit_transform(feat_train_scaled[scaling_features])
        
        # join training labels with scaled training features 
        labels_feat_train = pd.merge(labels_train, feat_train_scaled, on=id_string, how='inner')

        
    elif scaling == 'robust':
        # create RobustScaler
        rob_scal = RobustScaler(with_centering=True, with_scaling=True)

        # scale features
        feat_train_scaled = feat_train_filled.copy()
        feat_train_scaled[scaling_features] = rob_scal.fit_transform(feat_train_scaled[scaling_features])
        
        # join training labels with scaled training features 
        labels_feat_train = pd.merge(labels_train, feat_train_scaled, on=id_string, how='inner')
        
    elif scaling == None:
        labels_feat_train = pd.merge(labels_train, feat_train_filled, on=id_string, how='inner')
    
    else:
        raise ValueError("enter valid scaling ('standard', 'robust', None)")

    
    # perform PCA if number of pca components are given - TODO: omit?? => then delete pca_components argument in all run_training files
    if pca_components:
        if scaling != 'standard':
            raise ValueError("scaling should be 'standard' when performing PCA.")
        
        X = labels_feat_train.drop(columns=[id_string, label_string])
        pca = PCA(n_components=pca_components)   
        # get X_train and y_train as numpy arrays
        X_train = pca.fit_transform(X)
        y_train = labels_feat_train[label_string].values
        return X_train, y_train
        
    else:
        # get X_train and y_train as numpy arrays
        X_train = labels_feat_train.drop(columns=[id_string, label_string]).values
        y_train = labels_feat_train[label_string].values
        return X_train, y_train, features_selected
    

    
def get_feat_reduced_trainig_data(perc_var_exp,    
                                  model, 
                                  omit_scaling_features,
                                  scaling,
                                  score
                                 ):
    """
    Function to return training data on which feature reduction has been performed using PCA.
    
    """
    # get score specific arguments
    score_string, id_string, label_string = get_score_specific_args(score)
    
    # load reduced_training_data
    feat_train_red = pd.read_csv(f"{score_string}/training/train_test_data/feat_train_reduced_{perc_var_exp}_{config_vars[f'data_prep_date_{score}']}.csv.gz")
    labels_train = pd.read_csv(f"{score_string}/training/train_test_data/labels_train_{config_vars[f'data_prep_date_{score}']}.csv.gz")
    
    # fill features missing values
    feat_train_red_filled = fill_missing_vals(feat_train_red, model, score)
    
    # get features that should be scaled
    omit_scaling = get_features_from_groups(omit_scaling_features, feat_train_red, score) # features that should not be scaled
    scaling_features = [i for i in feat_train_red.columns if i not in omit_scaling + [id_string]] # features that should be scaled
    
    if scaling == 'standard':
        # create StandardScaler
        stand_scal = StandardScaler()

        # scale features
        feat_train_red_scaled = feat_train_red_filled.copy()
        feat_train_red_scaled[scaling_features] = stand_scal.fit_transform(feat_train_red_scaled[scaling_features])
        
        # join training labels with scaled training features 
        labels_feat_red_train = pd.merge(labels_train, feat_train_red_scaled, on=id_string, how='inner')

        
    elif scaling == 'robust':
        # create RobustScaler
        rob_scal = RobustScaler(with_centering=True, with_scaling=True)

        # scale features
        feat_train_red_scaled = feat_train_red_filled.copy()
        feat_train_red_scaled[scaling_features] = rob_scal.fit_transform(feat_train_red_scaled[scaling_features])
        
        # join training labels with scaled training features 
        labels_feat_red_train = pd.merge(labels_train, feat_train_red_scaled, on=id_string, how='inner')
        
    elif scaling == None:
        labels_feat_red_train = pd.merge(labels_train, feat_train_red_filled, on=id_string, how='inner')
    
    else:
        raise ValueError("enter valid scaling ('standard', 'robust', None)")

        
    # get X_train and y_train as numpy arrays
    X_train = labels_feat_red_train.drop(columns=[id_string, label_string]).values
    y_train = labels_feat_red_train[label_string].values
    
    features = [i for i in feat_train_red if i not in [id_string]]
    
    return X_train, y_train, features
           

    
# function to train using grid search and cross-validation
def train_with_grid_search(ID,
                           X_train, 
                           y_train, 
                           estimator, 
                           param_grid, 
                           cv, 
                           scoring, #='roc_auc', 
                           verbose,
                           smote,
                           #add_to_filename, #=''
                           score
                          ):
    
    # get score specific arguments
    score_string, _, _ = get_score_specific_args(score)
    
    # training starting time 
    start_time = datetime.today().strftime('%Y-%m-%d--%H-%M-%S')
    print(f"start_time: {start_time}")
    
    # create a GridSearchCV object with X-fold cross-validation
    grid_search = GridSearchCV(estimator=estimator, 
                               param_grid=param_grid, 
                               cv=cv, 
                               scoring=scoring, 
                               verbose=verbose)
    
    if smote:
        # add prefix to params in param grid (needed for pipeline)
        param_grid = {'classifier__' + key: value for key, value in param_grid.items()}

        # define SMOTE object
        smote_obj = SMOTE(sampling_strategy='auto', random_state=1)  # You can adjust parameters as needed

        # create a pipeline with SMOTE object and the classifier
        pipeline = Pipeline([('smote', smote_obj), ('classifier', estimator)])

        # create a GridSearchCV object with the imbalanced-learn pipeline
        grid_search = GridSearchCV(estimator=pipeline,
                                   param_grid=param_grid,
                                   cv=cv,  
                                   scoring=scoring,
                                   verbose=2)
 

    # fit the model to the training data
    grid_search.fit(X_train, y_train)
    
    # training ending time
    end_time = datetime.today().strftime('%Y-%m-%d--%H-%M-%S')
    print(f"end_time: {end_time}")

    # get the best parameters and estimator
    best_params = grid_search.best_params_
    best_classifier = grid_search.best_estimator_
    print('Best params: ', best_params)
    
    # get results df and add additional parameters
    cv_results = pd.DataFrame(grid_search.cv_results_)
    cv_results['start_time'] = start_time
    cv_results['end_time'] = end_time
    cv_results['ID'] = ID
    
    # write csv
    file_name = f"{score_string}/training/results/cv_results_ID{ID}.csv"
    cv_results.to_csv(file_name, index=False)
    
    
    results_dic = {
                  'cv_results': cv_results,
                  'best_params': best_params,
                  'best_classifier': best_classifier
                  }
    
    with open(f"{score_string}/training/results/results_ID{ID}.pkl", 'wb') as file:
        pickle.dump(results_dic, file)    

    return cv_results, best_params, best_classifier




def get_config_results_dics(ID, score):
    # get score specific arguments
    score_string, _, _ = get_score_specific_args(score)

    # load configuration dictionary
    with open(f"{score_string}/training/config_files/config_dic_ID{ID}.pkl", 'rb') as file:
        config_dic = pickle.load(file)
        
    # load results dictionary
    with open(f"{score_string}/training/results/results_ID{ID}.pkl", 'rb') as file:
        results_dic = pickle.load(file)
    
    return config_dic, results_dic


# function to get params in param_grid that have multiple values (for display in a heatmap)
def get_heatmap_params(param_grid):
    len_list = [len(param_grid[i]) for i in param_grid.keys()]
    indices = [i for i, element in enumerate(len_list) if element > 1]
    params = [list(param_grid.keys())[i] for i in indices]
    return params



# function to create a 2D heatmap with axes=hyperparamters and values=test score
def plot_2D_heatmap(ID, cv_results, param1, param2, figsize, save, show, score):
    
    # get score specific arguments
    score_string, _, _ = get_score_specific_args(score)
    
    # pivot the DataFrame to create a heatmap-ready format
    heatmap_data = cv_results.pivot(index='param_'+ param1 , 
                                    columns='param_'+ param2 , 
                                    values='mean_test_score')

    # create a figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    # create a heatmap
    cax = ax.matshow(heatmap_data, cmap='viridis')

    # set the tick labels for C and gamma
    ax.set_xticks(range(len(heatmap_data.columns)))
    ax.set_yticks(range(len(heatmap_data.index)))
    
    if is_numeric(heatmap_data.columns.values):
        ax.set_xticklabels(['{:.1e}'.format(i) for i in heatmap_data.columns.values])
    else:
        ax.set_xticklabels([i for i in heatmap_data.columns.values])
        
    if is_numeric(heatmap_data.index.values):
        ax.set_yticklabels(['{:.1e}'.format(i) for i in heatmap_data.index.values])
    else:
        ax.set_yticklabels([i for i in heatmap_data.index.values])

    # label the axes
    ax.set_xlabel(param2)
    ax.set_ylabel(param1)

    # add colorbar
    cbar = fig.colorbar(cax)

    # annotate each cell with score values
    for i in range(len(heatmap_data.index)):
        for j in range(len(heatmap_data.columns)):
            score = heatmap_data.values[i, j]
            ax.text(j, i, f'{score:.3f}', ha='center', va='center', color='w')

    # save figure
    if save:
        plt.savefig(f"{score_string}/training/results/heatmap2D_ID{ID}.png", format="png")
    
    
    # display the heatmap
    if show:
        plt.show()
 
    
    

# function to get feature importance by using permutation importance 
def get_permutation_importance(ID,
                               X_train, 
                               y_train, 
                               features,
                               classifier,
                               scoring,
                               plot,
                               random_state,
                               score
                              ):
    
    # get score specific arguments
    score_string, _, _ = get_score_specific_args(score)
    
    # fit best classfier with full training data
    classifier.fit(X_train, y_train)
    
    # get permutation importance
    perm_importance = permutation_importance(classifier, X_train, y_train, scoring=scoring, random_state=random_state)
    
    # create a df sorted by feature importance
    sorted_idx = perm_importance.importances_mean.argsort()
    feature_imp = pd.DataFrame({'feature': [features[i] for i in sorted_idx],
                                'perm_imp': perm_importance.importances_mean[sorted_idx]}).sort_values('perm_imp', ascending=False)
        
    feature_imp['scoring'] = scoring
    
    # save results
    date_time = datetime.today().strftime('%Y-%m-%d--%H-%M-%S')
    feature_imp.to_csv(f"{score_string}/training/results/perm_importance_ID{ID}_rs-{random_state}_{date_time}.csv", index=False)
    
    # plot feature importance
    if plot:
        fig = plt.figure(figsize=(6, 25))  
        plt.barh([features[i] for i in sorted_idx], perm_importance.importances_mean[sorted_idx])
        plt.xlabel("Permutation Importance")
        plt.show
    
    return(feature_imp)


def plot_2D_heatmap_fixed_params(ID, cv_results, best_params, param1, param2, figsize, save, show, score):
    """
    Function to plot a 2D heatmap of the given parameters param1 and param2. 
    The remaining parameters are fixed to their best values.
    """
    
    # get score specific arguments
    score_string, _, _ = get_score_specific_args(score)
    
    # get paramters that should be fixed to their best values
    param_cols = [i for i in cv_results.columns if i.startswith("param_")]
    fixed_params = {key: value for key, value in best_params.items() if key != param1 and key != param2}
    
    # get results dataframe with fixed parameters
    cv_results_fixed = cv_results.copy()
    for col, value in fixed_params.items():
        cv_results_fixed = cv_results_fixed[cv_results_fixed['param_' + col] == value]

    # pivot results dataframe
    heatmap_data = cv_results_fixed.pivot(index='param_'+ param1, 
                                    columns='param_'+ param2, 
                                    values='mean_test_score')

    # create a figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    # create a heatmap
    cax = ax.matshow(heatmap_data, cmap='viridis')

    # set the tick labels 
    ax.set_xticks(range(len(heatmap_data.columns)))
    ax.set_yticks(range(len(heatmap_data.index)))

    if is_numeric(heatmap_data.columns.values):
        ax.set_xticklabels(['{:.1e}'.format(i) for i in heatmap_data.columns.values])
    else:
        ax.set_xticklabels([i for i in heatmap_data.columns.values])

    if is_numeric(heatmap_data.index.values):
        ax.set_yticklabels(['{:.1e}'.format(i) for i in heatmap_data.index.values])
    else:
        ax.set_yticklabels([i for i in heatmap_data.index.values])

    # label the axes
    ax.set_xlabel(param2)
    ax.set_ylabel(param1)

    # add colorbar
    cbar = fig.colorbar(cax)

    # annotate each cell with score values
    for i in range(len(heatmap_data.index)):
        for j in range(len(heatmap_data.columns)):
            score = heatmap_data.values[i, j] * 100
            ax.text(j, i, f'{score:.1f}', ha='center', va='center', color='w')

    # add fixed parameters
    ax.text(0, len(heatmap_data.index) + 1, f'fixed parameters: {fixed_params}', ha='left', va='center', color='black')  
    
    # save figure
    if save:
        plt.savefig(f"{score_string}/training/results/heatmap2D_fixed_params_ID{ID}.png", format="png")
    
    # display the heatmap
    if show:
        plt.show()



def select_feat_from_perm_imp(ID, index, score):
    """
    Function to returns a list of the most important features according to permutation importance analysis
    (up to the given index) and a list of the remaining features.
    The index starts at 1.
    """
    # get score specific arguments
    score_string, _, _ = get_score_specific_args(score)
    
    # get permutation importance csv for given ID
    pattern = f"perm_imp_sum_rank_ID{ID}"
    files = os.listdir(f"{score_string}/training/results")
    matching_files = [file for file in files if file.startswith(pattern)]

    if len(matching_files) > 1:
        raise ValueError(f"Error: More than one permutation importance file for ID{ID}.")

    perm_imp = pd.read_csv(f"{score_string}/training/results/{matching_files[0]}").sort_values(by="rank_value", ascending=False)
    
    # select only the most important features up to given index
    important_feat = list(perm_imp.loc[0:index-1]['feature'])
    
    # get the remaining features
    unimportant_feat = list(perm_imp.loc[index:]['feature'])
    
    return important_feat, unimportant_feat



def png_to_pdf(png_files, output_pdf):
    # create a PDF file to save the images
    with PdfPages(output_pdf) as pdf:
        for png_file in png_files:
            # read and plot the PNG image using matplotlib
            img = plt.imread(png_file)

            # create a new figure and plot the image
            plt.figure(figsize=(img.shape[1] / 100, img.shape[0] / 100), dpi=100)
            plt.imshow(img)
            plt.axis('off')

            # save the figure to the PDF
            pdf.savefig()

            # close the figure
            plt.close()
            

def PCA_on_feature_group(group, 
                         model, 
                         n_components, 
                         plot,
                         omit_scaling_groups,
                        score):
    """
    Function to perform PCA on a specific feature group. 
    If plot = True: cumulative sum of variance is plotted.
    If n_components is an integer: a dataframe of the training data with the specified number of PC's is returned.
    """
    
    # get score specific arguments
    score_string, _, _ = get_score_specific_args(score)
    
    # load training features
    feat_train = pd.read_csv(f"{score_string}/training/train_test_data/feat_train_{config_vars[f'data_prep_date_{score}']}.csv.gz")
    
    # get features belonging to given feature group
    features_selected = get_features_from_groups([group], feat_train, score)
    
    # get training dataframe with only features from selected group
    feat_group_train = feat_train[features_selected]
    
    # fill missing values according to given model
    feat_group_train_filled = fill_missing_vals(df=feat_group_train, model=model, score=score)
    
    # create StandardScaler
    stand_scal = StandardScaler()
    
    # get features that should not be scaled
    omit_scaling = get_features_from_groups(['mgi', 'paralogues'], feat_train, score) 
    
    # get features that should be scaled
    scaling_features = [i for i in features_selected if i not in omit_scaling] 

    # scale features with StandardScaler
    feat_group_train_scaled = feat_group_train_filled.copy()
    feat_group_train_scaled[scaling_features] = stand_scal.fit_transform(feat_group_train_scaled[scaling_features])

    # perform PCA
    all_components = len(features_selected)

    X = feat_group_train_scaled

    full_pca = PCA(n_components=all_components)
    full_pca.fit(X)
    cumsum = np.cumsum(full_pca.explained_variance_ratio_)
    
    if plot:
        print(cumsum)
        plt.plot(cumsum, label = group)
        plt.xlabel("No. principal components")
        plt.ylabel("Explained Variance Ratio")
        plt.legend()

    if n_components:
        # get a dataframe with PC's
        pca = PCA(n_components=n_components)
        pca.fit(X)

        colnames = [f'{group}_PC{i+1}' for i in range(n_components)]
        pca_df = pd.DataFrame(pca.fit_transform(X), columns=colnames)
        pca_df['hgnc_id_int'] = feat_train['hgnc_id_int']
        
        # get loadings
        loadings = pca.components_
        loadings_df = pd.DataFrame(loadings, columns=X.columns, index=[f'PC{i+1}' for i in range(n_components)])
        
        return(pca_df, loadings_df, cumsum)
    else:
        return(None, None, cumsum)
    