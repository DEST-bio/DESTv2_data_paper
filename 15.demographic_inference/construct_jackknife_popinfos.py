import pandas as pd
import os

"""
Users will have to change `data_dir` to fit their machines. 
This script precedes `construct_jackknife_sfs.py`, which constructs site frequency
spectra from the popinfos written by this script. 
It is used to construct jackknifed replicates of popinfo files by randomly
removing one sample from each population with more than one sample.
"""

def main():
    # Setup
    data_dir = "data/"
    popinfo_dir = data_dir + "popinfos/"
    reps = 40 # number of trios to sample

    # Construct jackknife popinfos for all initial popinfos in popinfo_dir
    for popinfo_name in os.listdir(popinfo_dir):
        # Skip over non-popinfo files, like the folders created by this script
        if popinfo_name.endswith(".popinfo"):
            popinfo_name_no_ext = popinfo_name.rsplit(".", 1)[0]

            # Skip over popinfo files that already have jackknife replicates
            if os.path.exists(popinfo_dir + popinfo_name_no_ext + "_jackknife/"):
                print(f"Jackknife replicates of {popinfo_name_no_ext}.popinfo "
                      f"already exist. Skipping to next popinfo file.")
                continue

            # Print progress
            print(f"Constructing jackknife replicates of "
                  f"{popinfo_name_no_ext}.popinfo")

            construct_jackknife_popinfo(popinfo_dir, popinfo_name_no_ext, reps)

def construct_jackknife_popinfo(popinfo_dir: str, 
                                popinfo_name_no_ext: str, 
                                reps: int) -> None:
    """
    :param str popinfo_dir: Path to directory containing input popinfo file that
        that is to be jackknifed
    :param str popinfo_name_no_ext: Name of input input popinfo file, without its 
        extension
    :param int reps: Number of jackknife replicates to create
    :return: This function is void because it writes `reps` many jackknifed replicates
        of the file whose name, sans extension, is `popinfo_name_no_ext` and is
        located in `popinfo_dir`.
    :rtype: None
    """

    # Construct paths to input file and output directory
    input_popinfo = popinfo_dir + popinfo_name_no_ext + ".popinfo"
    output_dir = popinfo_dir + popinfo_name_no_ext + "_jackknife/"

    # Read in source popinfo
    popinfo = pd.read_csv(input_popinfo, sep="\t", header=None,
                          names=['sample', 'pop'])
    
    # Create output directory if it doesn't already exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Perform jackknifing replicates
    for rep in range(reps):
        # For pop-based groups with more than one member, exclude one random row
        jackknife_popinfo = popinfo.groupby('pop', group_keys=False).\
            apply(lambda group: group.sample(n=len(group)-1) if len(group) > 1 \
                                else group)
        
        # Write jackknifed popinfo to file
        jackknife_popinfo.to_csv(f"{output_dir}{popinfo_name_no_ext}_j{rep}.popinfo",
                                 sep="\t", header=False, index=False)

if __name__ == "__main__":
    main()
