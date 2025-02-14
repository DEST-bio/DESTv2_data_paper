## This directory contains David Bass's contributions to the project. It is still in rough shape and needs additional documentation and organization.

### Pipeline outline
- `convert_GDS_to_VCF.R`: Convert raw data in GDS file to VCF file
- `get_average_n_eff_per_samp.py`: Parse metadata and VCF to calculate average effective sample size (`n_eff`) per sample
- For each region:
  - `construct_popinfos.R`: Parse metadata to construct popinfo file describing each region (e.g. Europe, Transatlantic, Americas, mainland)
  - `get_average_n_eff_per_popinfo.py`: Parse output of `get_average_n_eff_per_samp.py` per popinfo
  - `construct_sfs.py`: Parse VCF and popinfo to construct SFSs
  - `run_moments.py`: Given an SFS and model, run moments to fit the model to the SFS as data
  - For jackknifing:
    - `construct_jackknife_popinfos.py`: Jackknife popinfo file, i.e. create repeats in which one locality is removed from each population with >1 sample
    - `construct_jackknife_sfs.py`: Parse VCF and jackknifed popinfo to construct jackknifed SFSs
- `jupyter_nbs/jackknife_inference.ipynb`: Perform analysis and visualization of `moments` demographic inference runs on jackknife-replicate SFSs.
- `get_Fsts.py`: Model selections with $F_st$, the classical population genetic statistic for measuring population differentiation.

### Demographic inference on DEST2 data.

![plot](figures/DEST2_fig_S16.png)

For split model, (est0, est1, est2, est3) is (N_EUW, N_EUE, T, m).

For admixture model, (est0, est1, est2, est3, est4, est5, est6, est7) is (N_1, N_2, N_admix, T_split, T_admix, m_2, m_3, p_admix), where N_1 describes the easternmost non-admixed population, N_2 describes the westernmost non-admixed population, m_2 is the migration rate in the 2-population epoch, and m_3 is the migration rate in the 3-population epoch.

To calculate collapsed log-likelihood, we sum 3D SFSs along their axis corresponding to SZ, “collapsing” them into a 2D SFS with the same shape as the SFSs used to fit the split model, then apply moments’ log-likelihood calculation function to collapsed data and model SFSs. The collapsed SFS effectively “integrates over” the alleles discovered in SZ, making it both geometrically and biologically comparable to the natively 2D SFSs yielded by the split model. This is a novel metric and is most appropriate for comparison of models with distinct numbers of populations.

The collpased likelihoods of all three admixture models far exceed the likelihood of the split model, indicating that 3-population models better describe 2-population models. It is worth noting that we have not searched the entire space of 2- and 3-population models, though doing so would likely be worthless because it is vast and composed of models that are biologically less intuitive than the admixture, split, and two_splits models that we have attempted to fit.

After rejecting the null hypothesis that a 2-population model is sufficient to describe European flies, we conclude that the admixture model in which SZ is the admixed population achieves the greatest likelihood. It is worth noting that collapsed likelihood favors the model in which EUW is the admixed population, though comparison via likelihoods, which are calculated on all available data, in the form of 3-dimensional SFSs, as opposed to down-projected, "collapsed" 2-dimensional SFSs.

Further analysis must investigate the potential of runaway behavior in certain parameter estimates, as it is implied to exist by migration rates of 0 and 10 (an upper bound provided to moments' optimization routine) in the EUW admixed and SZ admixed models, respectively.

We exclude two_splits models from this analysis because no convergence was observed in any of the three variants, indicating that, it is likely impossible to achieve a good fit to the data, and if a fit does exist, then it is likely rendered not biologically meaningful by extreme estimate values.

### Descriptions of each code file (roughly in order of use)
`convert_GDS_to_VCF.R` is an R script that converts the GDS file encoding all variants discovered in the DEST v2.0 dataset into a gzip'ed VCF file.

`get_average_n_eff_per_samp.py` is a Python script 

`get_average_n_eff_per_popinfo.py` is a Python script that takes in the output of `get_average_n_eff_per_samp.py`, stored in `metadata/mean_n_eff_per_sample_on_mainChroms.csv`, and the name (with extension) of a [popinfo](https://github.com/MomentsLD/moments/blob/main/moments/Misc.py#L576) file stored in `popinfos/` as a command-line argument, and outputs into `metadata/mean_n_eff_per_popinfo.csv`. The output values must then be copied into the dictionary `popinfo_metadata` at the top of `construct_sfs.py`. These values are average effective population sizes (across all samples to which a model is fit, as listed in a popinfo file that guides SFS construction) that are used in the rounding of pool-seq allele frequencies to integer counts.

`construct_popinfos.py`

`construct_sfs.py`

`wrapper_construct_sfs.slurm`

`construct_jackknife_popinfos.py`

`construct_jackknife_sfs.py`

`wrapper_construct_jackknife_sfs.slurm`

`moments_models.py` is a Python file containing functions that define the models used by `moments` for demographic inference. It is imported by `get_scaled_model_lls.py` and `run_moments.py`.

`run_moments.py`

`wrapper_run_moments.slurm`

`wrapper_run_moments_on_simulated_sfss.py`

`get_collected_output.py`

`jupyter_nbs/jackknife_inference_analysis.ipynb` is a Jupyter notebook containing final analysis and visualization of the output from `run_moments.py` on jackknife-replicate SFSs produced by `construct_jackknife_sfs.py`.

`get_optimal_estimate_confidence_ints` *this might be trash I need to double-check*

`get_Fsts.py`

`msprime_models.py` is a Python file containing functions that define the files used by `msprime` for demographic simulation. These models mirror the models in the `moments_models.py`, and they are used to perform simulation-based validation of collapsed-likelihood-based model selection procedure with `simulate_sfs.py` and `run_moments_on_simulated_sfss.py`. This simulation is reported in Text S3.

`simulate_sfs.py` is a Python script that uses `msprime`, together with demographic functions defined in `msprime_models.py`, to simulate replicates of American populations sampled in DEST v2.0 so that our collapsed-likelihood-based model selection procedure can be theoretically validated with `run_moments_on_simulated_sfss.py`. This simulation is reported in Text S3.

`wrapper_simulate_sfs.sh`

`run_moments_on_simulated_sfss.py` replicates the demographic inference performed on American populations in DEST v2.0, but on replicate data simulated with `msprime` in `simulate_sfs.py`. Its output is visualized and analyzed in `jupyter_nbs/CLL_validation.ipynb`.

`jupyter_nbs/CLL_validation.ipynb` visualizes and analyzes output from `run_moments_on_simulated_sfss.py` for Text S3.