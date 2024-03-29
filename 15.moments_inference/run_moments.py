import moments
import numpy as np
import sys
import moments_models as mm
import pandas as pd
import re

"""
This script is the heart of the analysis. It uses moments to fit a demographic model
to provided SFS and appends results to a central output file. It is meant to be
ran many times on each SFS in a parallel array job via the wrapper `wrapper_run_moments.slurm`.
"""

"""
Defining `pop_of_interest`: 
For the admixture model applied to Europe, `pop_of_interest` denotes the population
that is the result of the admixture event. 
For the twosplits model applied to Europe, `pop_of_interest` denotes the population
that is directly descended from the ancestral population, or, equivalently, the 
one that does not result from the splitting of an intermediate extinct population.
For the Transatlantic region, `pop_of_interest` values in the set {0, 1, 2, 3} respectively
encode the Americas being the result of admixture of Eastern Europe (EUE) and Guinea (GN),
Western Europe (EUW) and Guinea, EUE and Zambia (ZM), and EUW and ZM.
"""

# Global final variables
regions = ("Europe", "mainland", "Americas", "Transatlantic")
model_to_model_func = {"two_epoch": mm.two_epoch,
                       "split": mm.split,
                       "admixture": mm.admixture,
                       "twosplits": mm.twosplits}
model_to_param_names = {"two_epoch": ["nu", "T"],
                        "split": ["nu1", "nu2", "T", "m"],
                        "admixture": ["nu1", "nu2", "nu_admix", "T_split", 
                                      "T_admix", "m2", "m3", "admix_prop"],
                        "twosplits": ["nu1", "nu_intermediate", "nu2", "nu3", 
                                      "T1", "T2", "m2", "m3"]}
# Used to construct buffers of zeros in output such that all printed rows have the
# same num. of columns 
max_num_params = max([len(params) for params in model_to_param_names.values()])
# Map model name to "special population," which is the result of the
# admixture event in the `admixture` model and the population that is directly
# descended from the ancestral population in the `twosplits` model.
model_to_special_pop = {"admixture": 2, "twosplits": 0}

def main():
    # Handle command line arguments
    sfs_file = sys.argv[1]
    model = sys.argv[2]
    pop_of_interest = int(sys.argv[3]) if sys.argv[3].isnumeric() else "NA"
    output_file = sys.argv[4]
    maxiter = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
    opt_func_name = sys.argv[6] if len(sys.argv) > 6 else "Nelder-Mead"
    init_params_gen_mode = sys.argv[7] if len(sys.argv) > 7 else "uniform"
    
    # Setup
    region = get_region(sfs_file.split("/")[-1])

    # If the SFS is jackknifed, get its ID, which is between "_j" and ".npy"
    is_jackknife_sfs = bool(re.search("_j\d+.npy", sfs_file))
    jackknife_id = sfs_file.split("_j")[-1][:-4] if is_jackknife_sfs else "NA"

    run_moments(sfs_file, model, pop_of_interest, output_file, region, jackknife_id,
                maxiter, opt_func_name, init_params_gen_mode)

def run_moments(sfs_file: str,
                model: str,
                pop_of_interest,
                output_file: str,
                region: str,
                jackknife_id: str,
                maxiter: int, 
                opt_func_name: str, 
                init_params_gen_mode: str) -> None:
    """
    :param model str: Name of .npy file encoding SFS, with path
    :param model str: Name of demographic model to fit to SFS
    :param pop_of_interest: Integer index of population-of-interest in the SFS, 
        or "NA" if not applicable to model. The only models for which a population-
        of-interest is applicable are the three-population models, namely `admixture` 
        and `twosplits`.
    :param output_file str: Name of output file, including path
    :param region str: Region from which data was taken
    :param jackknife_id str: ID of jackknifing replicate used to choose which samples
        to include in the SFS
    :param maxiter int: Maximum number of iterations for which to run optimization
    :param opt_func_name str: Name of optimization function to use. Either "Nelder-Mead"
        or "BFGS", which are the two that moments implements wrappers for with the
        functions `optimize_log_fmin` and `optimize_log`, respectively.
    :param init_params_gen_mode str: Mode for generating initial parameters. Either

    TODO: Add last three last params to docstring 
    :return: None, writes line of moments output to file
    :rtype: None
    """

    model_func = model_to_model_func[model]

    if init_params_gen_mode == "from_collected_output":
        lower_bound, init_params, upper_bound = \
            get_moments_bounds_from_collected_output("output/collected_output.tsv", 
                                                    model, 
                                                    region, 
                                                    pop_of_interest)
    elif init_params_gen_mode == "uniform":
        lower_bound = {param: 1e-3 for param in model_to_param_names[model]}
        upper_bound = {param: val for param, val in zip(model_to_param_names[model],
                                                        [40, 40, 50, 50, 1, 0.1, 0.1, 1])}
        # `init_params` is intentionally not initialized here. This will be handled
        # in `run_optimization()`.
    else:
        raise ValueError(f"Initial parameter generation mode {init_params_gen_mode} "
                         "not recognized.")
    
    time_limit = 2 # in hours; this does nothing here, but is set in `wrapper_run_moments.slurm`
    
    if opt_func_name == "BFGS":
        opt_func = moments.Inference.optimize_log
    elif opt_func_name == "Nelder-Mead":
        opt_func = moments.Inference.optimize_log_fmin
    else:
        raise ValueError(f"Optimization function {opt_func_name} not recognized.")
    
    # Load SFS
    sfs = load_sfs(sfs_file, region, model, pop_of_interest)

    # Run moments
    init_params, opt_params, ll, func_calls, grad_calls = \
        run_optimization(sfs, lower_bound, upper_bound, model, model_func, 
                         opt_func=opt_func,
                         maxiter=maxiter,
                         verbose=0,
                         init_params=None)

    coll_pop_ll = get_collapsed_ll(sfs, opt_params, model_func, region, 
                                   pop_of_interest, model, ll)
        
    # Process output
    num_params = len(upper_bound)
    padding_zeros = [0] * (max_num_params - num_params)
    output_list = np.concatenate(([model, pop_of_interest], 
                                    init_params, padding_zeros,
                                    opt_params, padding_zeros,
                                    list(upper_bound.values()), padding_zeros,
                                    [ll, coll_pop_ll, 
                                    func_calls, grad_calls,
                                    maxiter, time_limit,
                                    jackknife_id, region]))

    # Write output
    write_output(output_list, output_file)

def get_moments_bounds_from_collected_output(collected_op_file: str, 
                                             model: str, 
                                             region: str, 
                                             pop_of_interest) -> tuple:
    """
    :param collected_op_file str: Path to file containing moments output data that
        has been cleaned and combined with `get_collected_output.py`
    :param model str: Name of demographic model fit to SFS
    :param region str: Region from which data was taken
    :param pop_of_interest: Code for population(s)-of-interest from the region in 
        to be represented in the SFS, or "NA" if not applicable to model
    """
    opt_params_df = get_opt_param_df(collected_op_file)
    pop_of_interest_suffix = "_" + str(pop_of_interest) if pop_of_interest != "NA" else ""
    num_params = len(model_to_param_names[model])
    opt_params = list(opt_params_df.\
                      loc[(region, model + pop_of_interest_suffix)])[2:2+num_params]
    init_params = {param: val for 
                   param, val in zip(model_to_param_names[model], opt_params)}
    # Double optimal parameters to get upper bounds
    upper_bound = {param: 2 * val for param, val in init_params.items()}
    # Make sure that admixture proportion is bounded above by 100%
    if "admix_prop" in upper_bound:
        upper_bound["admix_prop"] = 1

    lower_bound = {param: 1e-3 for param in upper_bound}

    return lower_bound, init_params, upper_bound

def get_opt_param_df(collected_op_file: str) -> pd.DataFrame:
    """
    :param collected_op_file str: Path to file containing moments output data that
        been cleaned and filtered to only contain maximum-likelihood runs for each
        region--model combination by `get_collected_output.py`
    :return: DataFrame containing moments output data with named columns
    :rtype: pd.DataFrame
    """
    return pd.read_csv(collected_op_file,
                       names=["ll", "coll_pop_ll"] + \
                             [f'est{i}' for i in range(max_num_params)],
                       index_col=[0, 1])

def get_region(sfs_file: str) -> str:
    """
    :param str sfs_file: Path to SFS file
    :return: Region from which SFS was generated, selected from the global final
        variable `regions`
    :rtype: str
    """
    # Regions is a global variable defined at the top of this file
    for region in regions:
        if region in sfs_file:
            return region
    raise "Region not found in SFS file name."

def load_sfs(sfs_file: str, 
             region: str,
             model: str, 
             pop_of_interest: int) -> moments.Spectrum:
    """
    
    """
    sfs = moments.Spectrum(np.load(sfs_file)).fold()

    # If appropriate, swap order of pops. to facilitate collapsing of SFS. This
    # condition should only be satisfied for Europe.
    if region == "Europe" and model in model_to_special_pop:
        print("Swapping according to rules for Europe")
        sfs = sfs.swapaxes(pop_of_interest, model_to_special_pop[model])
    # This conditional block was for the prior, unsupported definition of the Transatlantic
    # regions, in which only the Caribbean population was assumed to exist.
    # if region == "Transatlantic":
    #     print("Swapping according to rules for Transatlantic")
    #     sfs = sfs.swap_axes(2, 3) # Swap mainland with Guinea so that mainland is last
    #     sfs = sfs.marginalize([pop_of_interest]) # excluded_pop in {0, 1, 2}
    if region == "Transatlantic":
        # 0 = 00 = (EUE, Guinea)
        # 1 = 01 = (EUW, Guinea)
        # 2 = 10 = (EUE, Zambia)
        # 3 = 11 = (EUW, Zambia)
        # pop_ids=["Americas", "Europe_east", "Europe_west", "Guinea", "Zambia"]
        print("Reformatting SFS according to rules for Transatlantic v2.")
        # Remove unused European and African populations from SFS
        unused_pop_lists = [[2, 4], [1, 4], [2, 3], [1, 3]]
        sfs = sfs.marginalize(unused_pop_lists[pop_of_interest])
        # Swap Americas with the African population to make it the admixed one.
        # Thus the order of axes in the SFS, is African, European, Americas.
        sfs = sfs.swapaxes(0, 2)

    return sfs

def get_init_params(lower_bound: dict, 
                    upper_bound: dict, 
                    model=None) -> dict: # TODO: Fix lazy duck typing of `model`
    """
    :param lower_bound dict: Map of parameter names to lower bound values. For initial
        probing of parameter spaces, all lower bounds are set to 1e-3.
    :param upper_bound dict: Map of parameter names to upper bound values, the primary
        knob to turn when adjusting the size of the parameter space to explore at
        first
    :param model: Name of the demographic model specified by the the params
    :return: Dictionary of initial parameters, uniformly sampled from within the
        element-wise bounds provided by `lower_bound` and `upper_bound`
    :rtype: dict
    """

    # Uniformly sample from within the bounds
    init_params = {param: val 
                   for param, val 
                   in zip(lower_bound, 
                          np.random.uniform(list(lower_bound.values()), 
                                            list(upper_bound.values())))}

    # In admixture and twosplits models the earlier time (T_split and T1) must be
    # greater than the later time (in coalescent units).
    if model == "admixture":
        if init_params['T_split'] < init_params['T_admix']:
            init_params['T_split'], init_params['T_admix'] = \
                init_params['T_admix'], init_params['T_split']
    elif model == "twosplits":
        if init_params['T1'] < init_params['T2']:
            init_params['T1'], init_params['T2'] = \
                init_params['T2'], init_params['T1']
    
    return init_params

def write_output(output_list: list, output_file: str) -> None:
    """
    :param output_list list: List of output values to be written to the output file
    :param output_file str: Name of output file, including path
    :return: None, writes list of outputs as tab-separated to file
    :rtype: None
    """
    with open(output_file, "a") as f:
        output_string = ""
        for output in output_list:
            output_string += str(output) + "\t"
        output_string = output_string.strip()
        output_string += "\n"

        f.write(output_string)

def minor_perturb(params_dict: dict, factor: float=0.2) -> np.ndarray:
    """
    :param param_dict dict: Dictionary of parameter names to their values
    :param factor float: Factor by which to perturb each parameter
    :return: List of values in `params_dict` uniformly randomly perturbed by a proportion
        up to `factor`
    :rtype: np.ndarray
    """
    arr = np.array(list(params_dict.values()))
    perturbed_arr = arr + np.random.uniform(arr * -factor, arr * factor)
    for i, param in enumerate(perturbed_arr):
        perturbed_arr[i] = max(1e-3, param)
    
    # Ensure that admixture proportion does not exceed 100%
    if "admix_prop" in params_dict:
        p_idx = list(params_dict.keys()).index("admix_prop")
        perturbed_arr[p_idx] = min(1, perturbed_arr[p_idx])
    
    return perturbed_arr

# Get log-likelihood of collapsed SFS, i.e. the SFS with the suture zone pop.
# summed over, yielding 2D-SFSs whose likelihood is directly comparable to those
# achieved by 2-pop. models, e.g. split.
def get_collapsed_ll(sfs, opt_params, model_func, region, pop_of_interest, model, ll):
    # Get collapsed-population log-likelihood for 3-pop models.
    # No such calculated is needed for Transatlantic models and one-population models.
    if region == "Europe" and model in ["admixture", "twosplits"]:
        suture_pop_axis = 1
        if pop_of_interest == suture_pop_axis:
            coll_pop_axis = model_to_special_pop[model]
        else:
            coll_pop_axis = 1
    elif region in ["mainland", "Americas"] and model == "split":
        coll_pop_axis = 1
    else:
        return ll

    model_sfs = model_func(opt_params, sfs.sample_sizes)
    coll_sfs, coll_model_sfs = [arr.marginalize([coll_pop_axis]) 
                                for arr in [sfs, model_sfs]]
    
    return moments.Inference.ll_multinom(coll_model_sfs, 
                                         coll_sfs / coll_sfs.S())

def run_optimization(sfs, lower_bound, upper_bound, model, model_func, opt_func,
                     verbose=0, maxiter=1000, init_params=None):
    """
    
    """
    # Get initial "guess" parameters from uniform distribution if not provided
    if not init_params:
        init_params = get_init_params(lower_bound, upper_bound, model)
        init_params_list = list(init_params.values())
    # If initial parameters are are provided, tweak each of them by up to 20%
    elif init_params:
        init_params_list = minor_perturb(init_params)

    # Run optimization. This is the step that takes a while!
    op = opt_func(init_params_list, 
                  sfs, 
                  model_func,
                  lower_bound=list(lower_bound.values()), 
                  upper_bound=list(upper_bound.values()), 
                  multinom=True,
                  full_output=True, 
                  verbose=verbose,
                  maxiter=maxiter)

    # Parse output. Description of output found with 
    # help(scipy.optimize.fmin_bfgs)
    opt_params = op[0]
    func_calls, grad_calls = [0, 0] # dummy values because previous values of this 
    # have been incorrect due to variable reporting from NM and BFGS that I failed to
    # account for

    ll = moments.Inference.ll_multinom(model_func(opt_params, sfs.sample_sizes), 
                                       sfs / sfs.S())
    
    return init_params_list, opt_params, ll, func_calls, grad_calls

if __name__ == "__main__":
    main()
