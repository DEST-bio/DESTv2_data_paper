import numpy as np
import functools
import sys
import msprime
import msprime_models as msm
import string
import pandas as pd

model_name_to_data = {"one_pop": {'func': msm.get_one_population_dem,
                                         'num_params': 2,
                                         'num_pops': 1},
                      "split": {'func': msm.get_split_dem,
                                'num_params': 4,
                                'num_pops': 2},
                      "admixture": {'func': msm.get_admixture_dem,
                                    'num_params': 8,
                                    'num_pops': 3},
                      "twosplits": {'func': msm.get_two_splits_dem,
                                     'num_params': 8,
                                     'num_pops': 3}}

def main():
    # Handle command line arguments
    region_name = "mainland"
    model_name = "one_pop"

    # Get optimal parameter estimates for model of region from file
    est_file = "output/opt_est_confidence_intervals.tsv"
    est_df = pd.read_csv(est_file, sep="\t")
    est_df = est_df[(est_df["region"] == region_name) & (est_df["model"] == model_name)]
    est_ci_lowers = est_df.iloc[:, 2::2].values[0]
    est_ci_uppers = est_df.iloc[:, 3::2].values[0]

    dem_params = [np.random.uniform(low=low, high=high) for low, high in zip(est_ci_lowers, est_ci_uppers)]

    # Translate command line arguments into arguments for `simulate_sfs()`
    model_func, num_params, num_pops = extract_model_metadata(model_name)
    dem_params = dem_params[:num_params]

    output_dir = "data/"
    reps = 1

    simulate_sfs(region_name, model_name, model_func, dem_params, num_pops, output_dir, reps)

def extract_model_metadata(model_name):
    if "_" in model_name:
        if model_name.rsplit("_", 1)[1].isnumeric():
            model_name_base = model_name.rsplit("_", 1)[0]
            return model_name_to_data[model_name_base].values()
    return model_name_to_data[model_name].values()    

def simulate_sfs(region_name, model_name, model_func, dem_params, num_pops, output_dir, reps):
    # Ancestral population size, from Table 1 of Kapopoulou et al. 2020 (doi:10.1038/s41598-020-79720-1)
    n_anc = 17734#4
    # Lengths of 2L, 2R, 3L, 3R chrs. from Table 1 of Hoskins et al. 2015 (doi:10.1101/gr.185579.114)
    seq_l = 23.5e6 + 25.3e6 + 28.1e6 + 32.1e6
    # Average recombination rate and mutation rates from abstract of Wang et al. 2023 (doi:10.1101/gr.277383.122)
    recomb = 0
    mut = 3.265e-7

    samp_N = 4
    chrom_N = samp_N * 2 # x2 for diploidy
    sample_nodes = [[i for i in range(samp_N)], [i for i in range(samp_N, samp_N * 2)]]
    print(sample_nodes)
    dem = model_func(n_anc, dem_params)
    ns = {pop_id: samp_N for pop_id in string.ascii_uppercase[:num_pops]}

    # Get TreeSequence object
    tss = msprime.sim_ancestry(num_replicates=reps,
                               samples=ns, 
                               demography=dem,
                               sequence_length=seq_l,
                               recombination_rate=recomb)
    # Simulate mutations
    mutations_func = functools.partial(msprime.sim_mutations, rate=mut)
    mutations_reps = map(mutations_func, tss)
    # Get SFS from mutations added to TreeSequence object
    get_sfs_func = functools.partial(get_sfs, sample_nodes=sample_nodes)
    sfs_reps = map(get_sfs_func, mutations_reps)

    # Save SFSs to file
    for i, sfs in enumerate(sfs_reps):
        print(sfs)
        output_file = f"{output_dir}{region_name}_{model_name}_r{i}"
        np.save(output_file, sfs)

def get_sfs(mts, sample_nodes):
    spectrum = mts.allele_frequency_spectrum(sample_sets=sample_nodes, polarised=False)
    sfs = spectrum.data
    sfs[tuple(0 for _ in sfs.shape)] = 0
    sfs /= np.sum(sfs)
    return sfs

if __name__ == "__main__":
    main()
