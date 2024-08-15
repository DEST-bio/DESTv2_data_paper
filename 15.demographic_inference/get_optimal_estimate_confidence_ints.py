import pandas as pd
import numpy as np

"""
This script takes output from fine-tuned moments model fitting on jackknife SFS
replicates and calculates the 95% confidence intervals for the parameter estimates
of each region's optimal model, as determined in `jackknife_inference_analysis.ipynb`.
The data cleaning is adapted from `jackknife_inference_analysis.ipynb`.
"""

def main():
    data_dir = "output/"
    output_file = data_dir + "opt_est_confidence_intervals.tsv"

    # Setup; data loading and cleaning
    data_file = data_dir + "moments_output_jackknife.tsv"
    max_num_params = 8

    get_optimal_estimate_confidence_ints(data_file, output_file, max_num_params)

def get_optimal_estimate_confidence_ints(data_file, output_file, max_num_params):
    est_cols = ['est' + str(i) for i in range(max_num_params)]

    # Read in data
    df = pd.read_csv(data_file, sep='\t', header=None,
                        names=['model', 'pop_of_interest'] + \
                            ['init' + str(i) for i in range(max_num_params)] + \
                            est_cols + \
                            ['upper_bound' + str(i) for i in range(max_num_params)] + \
                            ['log_likelihood', 'collapsed_pop_ll', 'func_calls', 'grad_calls',
                                'maxiter', 'hour_limit',
                                'jackknife_id', 'region'])

    # Filter out regions for which no CLLs were used because all models had the
    # same dimensionalities
    df = df[(df['region'] == "mainland") | \
            (df['region'] == "Americas") | \
            (df['region'] == "Europe")]


    # Get the best fit for each jackknife replicate
    df = df.loc[df.groupby(['region', 'model', 'pop_of_interest', 'jackknife_id'],
                        dropna=False)\
                ['log_likelihood'].idxmax()]
    df = df[['region', 'model', 'pop_of_interest', 'jackknife_id', 'log_likelihood', 'collapsed_pop_ll'] + \
            est_cols]

    df = df.replace({'pop_of_interest': {0.0: '0', 
                                        1.0 : '1',
                                        2.0 : '2',
                                        3.0 : '3',
                                        np.nan: 'NA'},
                    'ci_bound': {0: 'lower', 1: 'upper'}})

    # Merge model and pop of interest, which do not need to be distinguished now that
    # we're only considering models within each region, and not calling general model
    # functions.
    pop_of_interest_suffices = ['_' + poi if poi != 'NA' else '' for poi in df['pop_of_interest']]
    df['model'] = df['model'] + pop_of_interest_suffices
    df.drop(columns=['pop_of_interest'], inplace=True)

    #####
    # Everything below this line is clean-up of old models and model names
    #####

    # Filter out attempted, but unused split-with-asymmetric-migration model.
    df = df[df['model'] != 'split_asymmig']
    # Rename models for clarity
    df[df['region'] == "Transatlantic"] = df.replace({'model': {'admixture_0': 'EUExWAfr',
                                                                'admixture_1': 'EUWxWAfr',
                                                                'admixture_2': 'EUExSAfr',
                                                                'admixture_3': 'EUWxSAfr'}})
    df[df['region'] == "Australia"] = df.replace({'model': {'admixture_0': 'EUWxWAfr',
                                                            'admixture_1': 'EUExWAfr',
                                                            'admixture_2': 'EUWxSAfr',
                                                            'admixture_3': 'EUExSAfr'}})
    df = df.replace({'model': {'two_epoch': 'one_pop'}})

    # region_to_optimal_model = {"mainland": "one_pop",
    #                            "Americas": "one_pop",
    #                            "Europe": "split",
    #                            "Transatlantic": "EUExSAfr",
    #                            "Australia": "EUWxSAfr"}
    ests_df_arr = []
    #for region, model in region_to_optimal_model.items():
    for region in df['region'].unique():
        for model in df[df['region'] == region]['model'].unique():
            row = [region, model]
            for i in range(max_num_params):
                ci = np.array(df[(df.region == region) & (df.model == model)]['est' + str(i)].quantile([0.025, 0.975]))
                row.extend(ci)
            ests_df_arr.append(row)

    ests_df = pd.DataFrame(ests_df_arr, 
                       columns=['region', 'model'] + [f'est{i}_{bound}' 
                                                      for i in range(max_num_params) 
                                                      for bound in ['lower', 'upper']])
    ests_df.to_csv(output_file,
                   sep='\t', index=False)

if __name__ == "__main__":
    main()