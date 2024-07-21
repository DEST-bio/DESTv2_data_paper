import numpy as np
import functools
import sys
import msprime
import msprime_models as msm
import string

model_name_to_data = {"one_population": {'func': msm.get_one_population_dem,
                                         'num_params': 2,
                                         'num_pops': 1},
                      "split": {'func': msm.get_split_dem,
                                'num_params': 4,
                                'num_pops': 2},
                      "admixture": {'func': msm.get_admixture_dem,
                                    'num_params': 8,
                                    'num_pops': 3},
                      "two_splits": {'func': msm.get_two_splits_dem,
                                     'num_params': 8,
                                     'num_pops': 3}}

def main():
    # Handle command line arguments
    model_name = sys.argv[1]
    samp_N = sys.argv[2]
    output_dir = sys.argv[3]

    # Translate command line arguments into arguments for `simulate_sfs()`
    model_func, num_params, num_pops = model_name_to_data[model_name].values()

    simulate_sfs(model_func, num_params, num_pops, samp_N, output_dir)

def simulate_sfs(model_func, num_params, num_pops, samp_N, output_dir):
    # msprime demographic and ancestry simulation parameters (unchanged params)

    # Ancestral population size, from Table 1 of Kapopoulou et al. 2020 (doi:10.1038/s41598-020-79720-1)
    n_anc = 177344
    ns = {pop_id: samp_N for pop_id in string.ascii_uppercase[:num_pops]}
    # Lengths of 2L, 2R, 3L, 3R chrs. from Table 1 of Hoskins et al. 2015 (doi:10.1101/gr.185579.114)
    seq_l = 23.5e6 + 25.3e6 + 28.1e6 + 32.1e6
    # Average recombination rate and mutation rates from abstract of Wang et al. 2023 (doi:10.1101/gr.277383.122)
    recomb = 2.75e-8
    mut = 3.265e-9 # mutation rate

    # need to make list of sample node ids to convert from ts to joint afs
    chrom_N = samp_N * 2
    sample_nodes = [list(range(i, i + chrom_N)) for i in range(0, chrom_N * num_pops, chrom_N)]

    # Handle command line arguments
    param_id = sys.argv[1]
    reps = int(sys.argv[2])
    # dem_params = [nu1, nu2, nu_admix, T_split, T_admix, p_admix]
    dem_params = [float(param) for param in sys.argv[3:3+num_params]]

    dem = model_func(n_anc, dem_params)

    tss = msprime.sim_ancestry(num_replicates=reps,
                               samples=ns, demography=dem,
                               sequence_length=seq_l,
                               recombination_rate=recomb)

    mutations_func = functools.partial(msprime.sim_mutations, rate=mut)
    mutations_reps = map(mutations_func, tss)

    get_sfs_func = functools.partial(get_sfs, sample_nodes=sample_nodes)
    sfs_reps = map(get_sfs_func, mutations_reps)

    for i, sfs in enumerate(sfs_reps):
        output_file = output_dir + "p" + param_id + "r" + str(i)
        np.save(output_file, sfs)

def get_sfs(mts, sample_nodes):
    spectrum = mts.allele_frequency_spectrum(sample_sets=sample_nodes, polarised=False)
    sfs = spectrum.data
    sfs[tuple(0 for _ in sfs.shape)] = 0
    sfs /= np.sum(sfs)
    return sfs

if __name__ == "__main__":
    main()
