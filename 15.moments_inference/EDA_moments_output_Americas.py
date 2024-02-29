import numpy as np
import pandas as pd
from tabulate import tabulate

"""
This script, alongside all others named `EDA_moments_output_*.py`, is an undocumented 
mess of redundant code that was used to explore the output of initial the moments runs.
Its output was used to determine the starting parameters for moments runs on jackknifed
SFSs.
"""

# Load in data (stored in local Ubuntu subsystem)
dir = "output/"
max_num_params = 4
est_cols = ['est' + str(i) for i in range(max_num_params)]
df = pd.read_csv(dir + 'moments_output_Americas.tsv', sep='\t', header=None,
                 names=['model', 'pop_of_interest'] + \
                       ['init' + str(i) for i in range(max_num_params)] + \
                       est_cols + \
                       ['upper_bound' + str(i) for i in range(max_num_params)] + \
                       ['ll', 'coll_pop_ll', 'func_calls', 'grad_calls',
                        'maxiter', 'hour_limit'])

df['model_and_pop'] = df.model + '_' + df.pop_of_interest.astype(str)
model_and_pops = df.model_and_pop.unique()

# Given a Pandas series 'col' and a tolerance proportion 'tol_prop', return a boolean
# vector that equals the condition that each element is within tol_prop percent
# of the maximum value in col.
def is_clustered_close_to_max(col, tol_prop):
    tol = abs(col.max()) * tol_prop
    tol_threshold = col.max() - tol
    return col > tol_threshold
# Filter out low-likelihood runs
opt_ll_filt = df.groupby('model_and_pop')['ll'].transform(is_clustered_close_to_max,
                                                          tol_prop=0.001)

print(opt_ll_filt[opt_ll_filt.isna()].index)
df = df.loc[opt_ll_filt]

# Number of runs per model with log-likelihoods within 10% of the maximum 
# discovered log-likelihood for that model
print(df.groupby('model_and_pop').count().model)

df_opt = df.groupby('model_and_pop').max('ll')[['ll', 'coll_pop_ll'] +\
                                                        est_cols]

print()
print(df_opt)

df_opt_common = df_opt.copy()
df_opt_common['region'] = "Americas"
df_opt_common = df_opt_common[['region'] + list(df_opt_common.columns[:-1])]
for i in range(4, 8):
    df_opt_common[f'est{i}'] = 0
df_opt_common.to_csv("output/collected_output.tsv", mode='a', header=None)

print()
print(df[(df.model_and_pop != "split_nan")].\
      sort_values(['model', 'pop_of_interest'])\
      [['model_and_pop'] + \
       est_cols + \
       ['ll', 'coll_pop_ll']])

df_tab = df_opt.reset_index().copy()
df_tab['num_conv'] = df.groupby('model_and_pop').count().model.tolist()

1/0

df_tab = df_tab[~df_tab.model_and_pop.str.contains("twosplits")]
df_tab['model'] = df_tab['model_and_pop'].replace(['split_nan', 'admixture_0.0', 
                                                   'admixture_1.0', 'admixture_2.0'],
                                                   ['split', 'EUW admixed',
                                                    'SZ admixed', 'EUE admixed'])
df_tab = df_tab.round(2)
df_tab = df_tab[['model', 'num_conv', 'll', 'coll_pop_ll'] + est_cols]
df_tab.columns = ['Model', 'Num. converged', 'Max. log-likelihood',
                  'Max. collapsed log-likelihood'] + est_cols

df_tab = df_tab.loc[[3, 0, 1, 2]]

df_tab['Max. collapsed log-likelihood'].loc[3] = ""
for i in range(4, 8):
    df_tab[f'est{i}'].loc[3] = ""


print(df_tab.to_markdown(index=False))
