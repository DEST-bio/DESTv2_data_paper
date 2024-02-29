import itertools as it

"""
This script outputs an options file for a Slurm array job to run via `wrapper_run_moments.slurm`.
"""

def main():
    sfs_dir = "/project/berglandlab/david_DEST2/sfss/"
    options_file = "options_run_moments_Europe.txt"
    n_reps = 9999 // 7
    
    models = ["admixture", "twosplits"]
    moments_output_file = "output/moments_output_Europe.tsv"

    unique_line_arrs = [[sfs_dir + "k4_Europe.npy", "split", "NA", moments_output_file]]
    for sfsf, model, poi, mof in it.product([sfs_dir + "k8_Europe.npy"],
                                           models,
                                           range(3),
                                           [moments_output_file]):
        unique_line_arrs.append([sfsf, model, str(poi), mof])
    
    for line in unique_line_arrs:
        print(line, "\t", len(line))

    with open(options_file, 'w') as f:
        for sfsf, model, poi, mof in unique_line_arrs:
            for _ in range(n_reps):
                f.write("\t".join([sfsf, model, poi, mof]) + "\n")
        

if __name__ == "__main__":
    main()
