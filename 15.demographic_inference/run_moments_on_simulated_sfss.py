import sys
import simulate_sfs as ss
import run_moments as rm
import moments
import numpy as np

"""
This script runs moments on msprime-simulated SFSs in order to validate the strategy
of demographic model selection based on collapsed log-likelihoods (CLLs).
"""

def main():
    # Handle command line arguments (identical to `simulate_sfs.py`)
    region = sys.argv[1]
    sim_model = sys.argv[2]
    dem_rep = sys.argv[3]

    sfs_dir = "data/simulated_sfss/"
    bisect_name = "_bisected" if sim_model == "one_pop" else ""
    sfs_file_no_ext = f"{sfs_dir}{region}_{sim_model}{bisect_name}_d{dem_rep}_s"

    output_file = "output/moments_output_CLL_validation.tsv"
    sim_reps = 40

    # Get optimal parameter estimates for model of region from file
    est_file = "output/opt_est_confidence_intervals.tsv"
    est_ci_lowers, est_ci_uppers = ss.get_param_CIs(est_file, region, sim_model)
    est_lowers, est_uppers = double_est_cis(est_ci_lowers, est_ci_uppers)

    for sim_rep in range(sim_reps):
        for model in ["one_pop", "split"]:
            print(f"Sim rep: {sim_rep}, Model: {model}")
            sfs_file = f"{sfs_file_no_ext}{sim_rep}.npy"
            sfs = moments.Spectrum(np.load(sfs_file)).fold()

            if model == "one_pop":
                sfs = sfs.marginalize([1])

            output_list = rm.run_moments(sfs=sfs,
                                        model=model, 
                                        pop_of_interest="NA", 
                                        region=region,
                                        jackknife_id=f"{dem_rep}_{sim_rep}",
                                        maxiter=1000, 
                                        opt_func_name="BFGS", 
                                        init_params_gen_mode="uniform",
                                        lower_bound_vals=est_lowers,
                                        upper_bound_vals=est_uppers)
            rm.write_output(output_list, output_file)

def double_est_cis(est_ci_lowers, est_ci_uppers):
    est_ci_means = [(low + high) / 2 for low, high in zip(est_ci_lowers, est_ci_uppers)]
    est_ci_lowers = [max(low - abs(low - mean), 1e-3) 
                     for mean, low in zip(est_ci_means, est_ci_lowers)]
    est_ci_uppers = [high + abs(high - mean)
                     for mean, high in zip(est_ci_means, est_ci_uppers)]
    return est_ci_lowers, est_ci_uppers

if __name__ == "__main__":
    main()
