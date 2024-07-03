import itertools as it

"""
This script outputs an options file for a Slurm array job to run via `wrapper_run_moments.slurm`.
"""

def main():
    sfs_dir = "/project/berglandlab/david_DEST2/sfss/"
    options_file = "options_run_moments_Europe_asymmig.tsv"
    n_reps = 2500
    model = "split_asymmig"
    output_file = "output/moments_output_Europe.tsv"

    
    with open(options_file, "a") as f:
        for _ in range(n_reps):
            f.write("\t".join([f"{sfs_dir}k4_Europe.npy", model, "NA", output_file]) + "\n")

if __name__ == "__main__":
    main()


# /project/berglandlab/david_DEST2/sfss/k4_Europe.npy       split    NA      output/moments_output_Europe.tsv