import numpy as np
import pandas as pd
import run_moments as rm

def main():
    output_dir = "output/"

    for region in rm.regions:
        data_file = f"{output_dir}moments_output_{region}.tsv"
        max_num_params = region_to_max_num_params(region)
        write_data(data_file, max_num_params, region)

def region_to_max_num_params(region: str) -> int:
    if region in ['Europe', 'Transatlantic']:
        return 8
    elif region in ['mainland', 'Americas']:
        return 4
    raise "Invalid region."

def create_model_and_pop(row: pd.core.series.Series) -> str:
    """
    :param row pd.core.series.Series: Pandas DataFrame row object containing moments 
        output data
    :return: String containing both the model and population-of-interest, allowing
        for easy filtering by their combination
    
    """
    if not pd.isna(row['pop_of_interest']):
        return f"{row['model']}_{int(row['pop_of_interest'])}"
    else:
        return row['model']

def write_data(data_file: str, max_num_params: int, region: str) -> None:
    """
    :param data_file str: File containing moments output data, including path
    :param max_num_params int: Maximum number of parameters of any model applied
        to the region
    :param region str: Region whose SFS was analyzed by moments to produce the data
    :return: None, writes rows containing 
    :rtype: None
    """
    df = pd.read_csv(data_file, sep='\t', header=None,
                     names=['model', 'pop_of_interest'] + \
                           ['init' + str(i) for i in range(max_num_params)] + \
                           ['est' + str(i) for i in range(max_num_params)] + \
                           ['upper_bound' + str(i) for i in range(max_num_params)] + \
                           ['ll', 'coll_pop_ll', 'func_calls', 'grad_calls',
                            'maxiter', 'hour_limit'])
    
    # Merge `model` and `pop_of_interest` into a single column `model_and_pop`
    # for ease of filtering
    df['model_and_pop'] = df.apply(create_model_and_pop, axis=1)
    # Get row with greatest log likelihood for each model and population combination
    df = df.loc[df.groupby('model_and_pop')['ll'].idxmax()]
    # Keep only columns of interest
    df = df[['model_and_pop', 'll', 'coll_pop_ll'] + \
            [f'est{i}' for i in range(max_num_params)]]
    df['region'] = region
    df = df[['region'] + list(df.columns[:-1])]

    # Write to file
    df.to_csv("output/collected_output.tsv", mode='a', header=None, index=False)

if __name__ == "__main__":
    main()
