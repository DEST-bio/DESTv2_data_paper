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
    region_name = sys.argv[1]
    model_name = sys.argv[2]
    if len(sys.argv) == 4:
        bisect_pop = sys.argv[3] == "True"
    else:
        bisect_pop = False

    # Translate command line arguments into arguments for `simulate_sfs()`
    model_func, num_params, num_pops = extract_model_metadata(model_name)

    # Get optimal parameter estimates for model of region from file
    est_file = "output/opt_est_confidence_intervals.tsv"
    est_ci_lowers, est_ci_uppers = get_param_CIs(est_file, region_name, model_name)

    # Setup
    output_dir = "data/simulated_sfss/"
    dem_reps = 40
    sim_reps = 40

    for dem_rep in range(dem_reps):
        print(dem_rep)
        # Sample parameters from confidence intervals
        dem_params = [np.random.uniform(low=low, high=high) for low, high in zip(est_ci_lowers, est_ci_uppers)]
        dem_params = dem_params[:num_params]

        bisect_name = "_bisected" if bisect_pop else ""
        output_file_no_ext = f"{output_dir}{region_name}_{model_name}{bisect_name}_d{dem_rep}_s"

        simulate_sfs(model_func, dem_params, num_pops,
                     output_file_no_ext, sim_reps, bisect_pop)

def get_param_CIs(est_file, region_name, model_name):
    est_df = pd.read_csv(est_file, sep="\t")
    est_df = est_df[(est_df["region"] == region_name) & (est_df["model"] == model_name)]
    est_ci_lowers = est_df.iloc[:, 2::2].values[0]
    est_ci_uppers = est_df.iloc[:, 3::2].values[0]

    return est_ci_lowers, est_ci_uppers

def extract_model_metadata(model_name):
    if "_" in model_name:
        if model_name.rsplit("_", 1)[1].isnumeric():
            model_name_base = model_name.rsplit("_", 1)[0]
            return model_name_to_data[model_name_base].values()
    return model_name_to_data[model_name].values()

def simulate_sfs(model_func, dem_params, num_pops,
                 output_file_no_ext, sim_reps, bisect_pop):
    # Ancestral population size, from Table 1 of Kapopoulou et al. 2020 (doi:10.1038/s41598-020-79720-1)
    n_anc = 177344
    # Lengths of 2L, 2R, 3L, 3R chrs. from Table 1 of Hoskins et al. 2015 (doi:10.1101/gr.185579.114)
    seq_l = 23.5e6 + 25.3e6 + 28.1e6 + 32.1e6
    # Average recombination rate and mutation rates from abstract of Wang et al. 2023 (doi:10.1101/gr.277383.122)
    recomb = 2.75e-11 # turned down 3 orders of magnitude from 2.75e-8
    mut = 3.265e-9

    samp_N = 5
    chrom_N = samp_N * 2 # x2 for diploidy
    total_N = chrom_N * num_pops
    sample_nodes = [list(range(i, i + chrom_N))
                    for i in range(0, total_N * (bisect_pop + 1), chrom_N)]
    dem = model_func(n_anc, dem_params)
    ns = {pop_id: samp_N * (bisect_pop + 1) for pop_id in string.ascii_uppercase[:num_pops]}

    # Get TreeSequence object
    tss = msprime.sim_ancestry(num_replicates=sim_reps,
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
    for sim_rep, sfs in enumerate(sfs_reps):
        print(sim_rep, end=" ")
        np.save(f"{output_file_no_ext}{sim_rep}", sfs)
    print()

def get_sfs(mts, sample_nodes):
    spectrum = mts.allele_frequency_spectrum(sample_sets=sample_nodes, polarised=False)
    sfs = spectrum.data
    sfs[tuple(0 for _ in sfs.shape)] = 0
    sfs /= np.sum(sfs)
    return sfs

if __name__ == "__main__":
    main()
