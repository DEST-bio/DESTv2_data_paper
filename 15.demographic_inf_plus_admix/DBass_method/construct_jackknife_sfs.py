import sys
import construct_sfs as cs

"""
Users will have to change variables listed under the comment 'Setup' in `main()`
to fit the paths on their machines.
This script is run after `construct_jackknife_popinfos.py`, which writes the jackknifed
popinfos that this script uses to build site frequency spectra.
This script was run as a parallel array job via the wrapper `wrapper_construct_jackknife_sfs.slurm`. 
It is used to construct jackknifed replicates of SFS files by randomly
"""

def main():
    # Handle command line arguments
    popinfo_name_no_ext = sys.argv[1]
    region = popinfo_name_no_ext.rsplit("_", 1)[0] # Remove '_j{rep}' from end of name

    # Setup
    data_dir = "/scratch/djb3ve/DESTv2_data_paper/15.demographic_inference/data/"
    popinfo_dir = f"{data_dir}popinfos/{region}_jackknife/"
    sfs_dir = f"{data_dir}sfss/{region}_jackknife/"
    vcf_file = f"{data_dir}vcfs/dest2_clustered_mainChroms.vcf.gz"
    n = 20 # Number of samples to project down to

    construct_jackknife_sfs(popinfo_name_no_ext, region, popinfo_dir, sfs_dir, vcf_file, n)

def construct_jackknife_sfs(popinfo_name_no_ext: str, 
                            region: str,
                            popinfo_dir: str, 
                            sfs_dir: str, 
                            vcf_file: str, 
                            n: int) -> None:
    """

    """
    cs.make_dirs([popinfo_dir, sfs_dir])
    mean_n_eff, pop_ids = cs.popinfo_metadata[region]
    Npop = len(pop_ids)
    ns = [n] * Npop
    
    cs.write_sfs(vcf_file=vcf_file, 
                 popinfo_file=f"{popinfo_dir}{popinfo_name_no_ext}.popinfo",
                 pop_ids=pop_ids,
                 projections=ns,
                 n_eff=mean_n_eff,
                 rounding_method="counts_no_round_up",
                 sfs_file=f"{sfs_dir}{popinfo_name_no_ext}.npy")


if __name__ == "__main__":
    main()
    
