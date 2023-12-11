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

from itertools import combinations


# In[9]:


print(get_latest_ML_config_training_file(config_dir))
pd.read_csv(get_latest_ML_config_training_file(config_dir))


# In[3]:


ID = 37
config_dic, results_dic = get_config_results_dics(ID=ID)    


# In[4]:


print(results_dic['best_params'])
results_dic['cv_results']['mean_test_score'].max()


# In[7]:


# plot 2D heatmap with fixed other hyperparamters (to best values)    
plot_2D_heatmap_fixed_params(ID=ID, 
                            cv_results=results_dic['cv_results'], 
                            best_params = results_dic['best_params'],
                            param1='learning_rate', 
                            param2='n_estimators', 
                            figsize=(15,10), 
                            save=False, 
                            show=True)


# In[ ]:





# In[122]:


feature_groups = ['gnomad', 'cellxgene', 'descartes', 'gtex', 'mgi', 'paralogues', 'phasCons', 'CpG_o2e']


# In[ ]:





# In[33]:


ID = 17

pattern = f"perm_imp_sum_rank_ID{ID}"

files = os.listdir(results_dir)
matching_files = [file for file in files if file.startswith(pattern)]

if len(matching_files) > 1:
    raise ValueError(f"More than one permutation importance file for ID{ID}.")

# matching_files[0]

perm_imp = pd.read_csv(f"{results_dir}/{matching_files[0]}")
fig = plt.figure(figsize=(6, 25))  
plt.barh(perm_imp[::-1]['feature'], perm_imp[::-1]['rank_value'])
plt.xlabel("Permutation Importance")
plt.show


# In[37]:


list(perm_imp.loc[0:20]['feature'])


# In[18]:


pattern


# In[27]:


# feature importance
perm_imp = pd.read_csv("results/perm_imp_sum_rank_ID37_2023-11-14--14-04-34.csv")
fig = plt.figure(figsize=(6, 25))  
plt.barh(perm_imp[::-1]['feature'], perm_imp[::-1]['rank_value'])
plt.xlabel("Permutation Importance")
plt.show


# In[48]:


def select_feat_from_perm_imp(ID, index):
    """
    Function to returns a list of the most important features according to permutation importance analysis
    (up to the given index) and a list of the remaining features.
    The index starts at 1.
    """
    # get permutation importance csv for given ID
    pattern = f"perm_imp_sum_rank_ID{ID}"
    files = os.listdir(results_dir)
    matching_files = [file for file in files if file.startswith(pattern)]

    if len(matching_files) > 1:
        raise ValueError(f"Error: More than one permutation importance file for ID{ID}.")

    perm_imp = pd.read_csv(f"{results_dir}/{matching_files[0]}").sort_values(by="rank_value", ascending=False)
    
    # select only the most important features up to given index
    important_feat = list(perm_imp.loc[0:index-1]['feature'])
    
    # get the remaining features
    unimportant_feat = list(perm_imp.loc[index:]['feature'])
    
    return important_feat, unimportant_feat
    
    
important_feat, unimportant_feat = select_feat_from_perm_imp(ID=17, index=10)    
    


# In[50]:


unimportant_feat


# In[45]:


list(perm_imp.loc[0:5]['feature'])

perm_imp.loc[5:]['feature']


# In[41]:


def example_function(value):
    if value < 0:
        raise ValueError("Value should be a positive number")
    
    print("Value is:", value)
        
        
example_function(0)        


# In[ ]:





# In[23]:


# Correlation heatmap

group = ['gnomad']

sel_feat = get_features_from_groups(group, feat_train)
feat_train = pd.read_csv(f'{train_test_data_dir}/feat_train_{data_prep_date}.csv')
corr_matrix = feat_train[sel_feat].corr(method='spearman')

# create a heatmap
plt.figure(figsize=(15, 15))
heatmap = plt.matshow(corr_matrix, cmap='coolwarm')

# add colorbar
plt.colorbar(heatmap)

# set labels
plt.xticks(range(len(corr_matrix.columns)), corr_matrix.columns, rotation=90)
plt.yticks(range(len(corr_matrix.columns)), corr_matrix.columns)

# display the heatmap
plt.title('Correlation Heatmap')
plt.show



# In[270]:


data = feat_train[get_features_from_groups(group, feat_train)]

# generate all possible combinations of column pairs
column_pairs = list(combinations(data.columns, 2))

# calculate Spearman correlation for each pair and store results in a new DataFrame
correlation_data = []
for pair in column_pairs:
    col1, col2 = pair
    spearman_corr = data[[col1, col2]].corr(method='spearman').iloc[0, 1]
    correlation_data.append([col1, col2, spearman_corr])

# create a new DataFrame with correlation values
correlation_df = pd.DataFrame(correlation_data, columns=['col1', 'col2', 'corr_spearman'])
correlation_df['corr_spearman_abs'] = np.abs(correlation_df['corr_spearman'])
correlation_df = correlation_df.sort_values('corr_spearman_abs', ascending=False)


a = get_features_from_groups([group[0]], feat_train)
b = get_features_from_groups([group[1]], feat_train)

correlation_df.query("col1 in @a & col2 in @b", engine='python')



# In[200]:


b


# In[197]:


get_features_from_groups([group[0]], feat_train)


# In[136]:


correlation_df.columns


# In[63]:


# cv_results = pd.read_csv(f"results/cv_results_ID{ID}.csv")

# with open(f"{config_dir}/config_dic_ID{ID}.pkl", 'rb') as file:
#     config_dic = pickle.load(file)


# In[64]:


results_dic['cv_results']['mean_test_score'].max()


# In[107]:





# In[65]:


# results_dic['cv_results']


# In[ ]:





# In[138]:


with open(f"{results_dir}/results_ID{ID}.pkl", 'rb') as file:
    results_dic = pickle.load(file)
    
cv_results = results_dic['cv_results']
best_params = results_dic['best_params']

max_score = cv_results['mean_test_score'].max()
cv_results.query('mean_test_score == @max_score')


# In[247]:


max_score = cv_results['mean_test_score'].max()
cv_results.query('mean_test_score == @max_score')


# In[38]:


ID=15
cv_results = pd.read_csv(f"results/cv_results_ID{ID}.csv")

with open(f"{config_dir}/config_dic_ID{ID}.pkl", 'rb') as file:
    config_dic = pickle.load(file)
    
    
plot_2D_heatmap(ID=config_dic['ID'],
                cv_results=cv_results,
                param1=get_heatmap_params(param_grid=config_dic['param_grid'])[0],
                param2=get_heatmap_params(param_grid=config_dic['param_grid'])[1],
                figsize=( 15, 
                         15),
                save=False,
                show=True
               )
    
# config_dic    


# In[26]:


ID=1
cv_results = pd.read_csv(f"results/cv_results_ID{ID}.csv")

with open(f"{config_dir}/config_dic_ID{ID}.pkl", 'rb') as file:
    config_dic = pickle.load(file)
    
    
plot_2D_heatmap(ID=config_dic['ID'],
                cv_results=cv_results,
                param1=get_heatmap_params(param_grid=config_dic['param_grid'])[0],
                param2=get_heatmap_params(param_grid=config_dic['param_grid'])[1],
                figsize=( 15, 
                         15),
                save=False,
                show=True
               )


# In[299]:


def plot_2D_heatmap_fixed_params(ID, cv_results, best_params, param1, param2, figsize, save, show):
    """
    Function to plot a 2D heatmap of the given parameters param1 and param2. 
    The remaining parameters are fixed to their best values.
    """
    
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
        plt.savefig(f"{results_dir}/heatmap2D_fixed_params_ID{ID}.png", format="png")
    
    # display the heatmap
    if show:
        plt.show()
     
    
    
# plot_2D_heatmap_fixed_params(ID=16, 
#                             cv_results=cv_results, 
#                             param1='min_samples_leaf', 
#                             param2='min_samples_split', 
#                             figsize=(10,10), 
#                             save=False, 
#                             show=True)


# In[26]:


with open(f"{config_dir}/config_dic_ID{ID}.pkl", 'rb') as file:
    config_dic = pickle.load(file)


# In[8]:


# plot_2D_heatmap(ID=config_dic['ID'],
#                 cv_results=cv_results,
#                 param1=get_heatmap_params(param_grid=config_dic['param_grid'])[0],
#                 param2=get_heatmap_params(param_grid=config_dic['param_grid'])[1],
#                 figsize=( 15, 
#                          15),
#                 save=False,
#                 show=True
#                )


# In[45]:


pd.read_csv("results/perm_imp_sum_rank_ID23_2023-11-13--16-55-03.csv")


# In[43]:


ID = 23

config_dic, results_dic = get_config_results_dics(ID=ID)  

df_combined = pd.DataFrame()

for i in np.arange(5):
    print(i)
    perm_imp = get_permutation_importance(ID=ID,
                                          X_train=config_dic['X_train'], 
                                          y_train=config_dic['y_train'], 
                                          features=config_dic['features'], 
                                          classifier=results_dic['best_classifier'],
                                          scoring=config_dic['scoring'],
                                          plot=False,
                                          random_state=i
                                         )
    
    
    perm_imp['rank_value'] = 1 / perm_imp['perm_imp'].rank(ascending=False)
    
    df_combined = pd.concat([df_combined, perm_imp], ignore_index=True)
    print(len(df_combined))

sum_rank_values = df_combined.groupby('feature')['rank_value'].sum().reset_index().sort_values(by='rank_value', ascending=False)
    
sum_rank_values


# In[ ]:


ID = 23

config_dic, results_dic = get_config_results_dics(ID=ID)    


pi = get_permutation_importance(ID=ID,
                               X_train=config_dic['X_train'], 
                               y_train=config_dic['y_train'], 
                               features=config_dic['features'], 
                               classifier=results_dic['best_classifier'],
                               scoring=config_dic['scoring'],
                               plot=True,
                                random_state = 1
                              )










# In[20]:


pi['rank_value'] = 1 / pi['perm_imp'].rank(ascending=False)


# In[15]:


j = pi.copy()


# In[16]:


j['rank_value'] = 1 / j['perm_imp'].rank(ascending=False)


# In[17]:


j


# In[30]:





# In[7]:


# import pandas as pd
# import gseapy as gp

# # Example: Assuming you have pi1, pi2, pi3 DataFrames
# # Replace this with your actual DataFrames
# pi = pd.DataFrame({'feature': ['A', 'B', 'C'], 'perm_imp': [0.2, 0.4, 0.1]})
# pi2 = pd.DataFrame({'feature': ['B', 'D', 'E'], 'perm_imp': [0.5, 0.2, 0.3]})
# pi3 = pd.DataFrame({'feature': ['A', 'C', 'D'], 'perm_imp': [0.3, 0.1, 0.4]})

# # Combine the ranked lists into a single DataFrame
# # combined_ranked_list = pd.concat([pi[['feature, perm_imp']], pi2[['feature, perm_imp']], pi3[['feature, perm_imp']]], ignore_index=True)
# combined_ranked_list = pd.concat([pi, pi2, pi3], ignore_index=True)

# combined_ranked_list['perm_imp'] = combined_ranked_list['perm_imp'].astype(float)

# # Use gseapy for gene set enrichment analysis
# gene_sets = {'rankedList': combined_ranked_list['perm_imp'].values}

# # Perform GSEA
# results = gp.prerank(rnk=combined_ranked_list, gene_sets=gene_sets, min_size=1, max_size=len(combined_ranked_list))

# # Display the top enriched items
# top_enriched_items = results.res2d.sort_values('NES', ascending=False).head(10)
# print("Top Enriched Items:")
# print(top_enriched_items[['Term', 'NES', 'pval', 'fdr']])


# In[ ]:





# In[ ]:





# In[ ]:





# In[11]:


ID = 29

config_dic, results_dic = get_config_results_dics(ID=ID)    


pi = get_permutation_importance(ID=ID,
                               X_train=config_dic['X_train'], 
                               y_train=config_dic['y_train'], 
                               features=config_dic['features'], 
                               classifier=results_dic['best_classifier'],
                               scoring=config_dic['scoring'],
                               plot=True,
                                random_state = 1
                              )


# In[13]:


fig = plt.figure(figsize=(6, 25))  
plt.barh(pi[::-1]['feature'], pi[::-1]['perm_imp'])
plt.xlabel("Permutation Importance")
plt.show

