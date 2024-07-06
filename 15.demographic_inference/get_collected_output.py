import pandas as pd

region_to_max_num_params = {'Europe': 8, 
                            'Transatlantic_expandedAfr': 8, 
                            'Australia': 8,
                            'mainland': 4, 
                            'Americas': 4}

def main():
    output_dir = "output/"

    for region in region_to_max_num_params.keys():
        print(region)
        data_file = f"{output_dir}moments_output_{region}.tsv"
        max_num_params = region_to_max_num_params[region]
        write_data(data_file, max_num_params, region)

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
    names = ['model', 'pop_of_interest'] + \
            ['init' + str(i) for i in range(max_num_params)] + \
            ['est' + str(i) for i in range(max_num_params)] + \
            ['upper_bound' + str(i) for i in range(max_num_params)] + \
            ['ll', 'coll_pop_ll', 'func_calls', 'grad_calls',
            'maxiter', 'hour_limit', 'jackknife_id', 'region']

    df = pd.read_csv(data_file, sep='\s+', header=None,
                     names=names)
    
    
    # Merge `model` and `pop_of_interest` into a single column `model_and_pop`
    # for ease of filtering
    df['model_and_pop'] = df.apply(create_model_and_pop, axis=1)
    # Get row with greatest log likelihood for each model and population combination
    asdf = df.groupby('model_and_pop')['ll'].idxmax()
    print(2.11)

    print("--")
    print(df.head)
    print(df.shape)
    print(asdf.shape)
    print("--")

    print(asdf)

    df = df.loc[asdf]
    print(2.2)
    # Keep only columns of interest
    df = df[['model_and_pop', 'll', 'coll_pop_ll'] + \
            [f'est{i}' for i in range(max_num_params)]]
    df['region'] = region
    df = df[['region'] + list(df.columns[:-1])]

    # Write to file
    df.to_csv("output/collected_output.tsv", mode='a', header=None, index=False)

    print(4)

if __name__ == "__main__":
    main()
