import itertools as it

"""
This script outputs an options file for a Slurm array job to run via `wrapper_run_moments.slurm`
that runs moments on jackknife replicates.
"""

def main():
    # Setup
    data_dir = "data/sfss/"
    options_file = "options_run_moments_jackknife.tsv"
    jk_reps = 40
    n_reps = 16
    maxiter = 1000
    opt_func_name = "BFGS"
    init_params_gen_mode = "from_collected_output"
    
    moments_output_file = "output/moments_output_jackknife.tsv"

    # sfs_dirs = ["k4_Europe_jackknife/", "k8_Europe_jackknife/",
    #             "k8_Americas_p1_jackknife/", "k8_Americas_p2_jackknife/",
    #             "k8_mainland_p1_jackknife/", "k8_mainland_p2_jackknife/",
    #             "k8_Transatlantic_jackknife/"]
    # models_and_pops = [[["split", "NA"]],
    #                    [["admixture", "0"], ["admixture", "1"], ["admixture", "2"],
    #                     ["twosplits", "0"], ["twosplits", "1"], ["twosplits", "2"]],
    #                    [["two_epoch", "NA"]],
    #                    [["split", "NA"]],
    #                    [["two_epoch", "NA"]],
    #                    [["split", "NA"]],
    #                    [["admixture", "0"], ["admixture", "1"], ["admixture", "2"]]]

    # sfs_dirs = ["k4_Transatlantic_jackknife/"]
    # models_and_pops = [[["admixture", "0"], ["admixture", "1"], ["admixture", "2"], ["admixture", "3"]]]

    sfs_dirs = ["k4_Transatlantic_expandedAfr_jackknife/", 
                "k4_Australia_jackknife/",
                "k4_Europe_jackknife/"]
    models_and_pops = [[["admixture", "0"], ["admixture", "1"], ["admixture", "2"], ["admixture", "3"]],
                       [["admixture", "0"], ["admixture", "1"], ["admixture", "2"], ["admixture", "3"]],
                       [["split_asymmig", "NA"]]]
    

    with open(options_file, 'w') as f:
        for sfs_dir, model_and_pop_list in zip(sfs_dirs, models_and_pops):
            for model, pop in model_and_pop_list:
                for jk_rep in range(0, jk_reps):
                    for _ in range(n_reps):
                        sfs_file = f"{sfs_dir}{sfs_dir.rsplit('_', 1)[0]}_j{jk_rep}.npy"
                        f.write(f"{data_dir}{sfs_file}\t{model}\t{pop}\t{moments_output_file}\t"
                                f"{maxiter}\t{opt_func_name}\t{init_params_gen_mode}\n")

if __name__ == "__main__":
    main()