
"""
This script outputs an options file for a Slurm array job to run via `wrapper_run_moments.slurm`.
"""

def main():
    options_file = "options_run_moments_Transatlantic.tsv"
    model = "admixture"
    pop_of_interest_codes = range(4)
    moments_output_file = "output/moments_output_Transatlantic.tsv"
    data_dir = "data/sfss/"
    sfs_file = f"{data_dir}k4_Transatlantic.npy"
    reps = 9999 // len(pop_of_interest_codes)

    with open(options_file, "w") as f:
        for pop_of_interest in pop_of_interest_codes:
            for _ in range(reps):
                f.write("\t".join([sfs_file, model, str(pop_of_interest), 
                                   moments_output_file]) + "\n")

if __name__ == "__main__":
    main()
