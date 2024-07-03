
"""
This script outputs an options file for a Slurm array job to run via `wrapper_run_moments.slurm`.
"""

def main():
    data_dir = "/project/berglandlab/david_DEST2/sfss/"
    options_file = "options_run_moments_Americas.txt"
    n_reps = 9999 // (2 * 2) # x2 in denom bc I don't think that I need to do so many reps
    
    sfs_files = [f"{data_dir}k8_Americas_p{p}.npy" for p in range(1, 2+1)]
    models = ["two_epoch", "split"]
    moments_output_file = "output/moments_output_Americas.tsv"
    pop_of_interest = "NA"

    with open(options_file, 'w') as f:
        for sfs_file, model in zip(sfs_files, models):
            for _ in range(n_reps):
                f.write("\t".join([sfs_file, model, pop_of_interest, moments_output_file]) + "\n")
        

if __name__ == "__main__":
    main()
