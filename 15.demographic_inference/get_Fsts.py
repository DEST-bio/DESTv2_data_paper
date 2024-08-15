import moments
import numpy as np

def main():
    data_dir = "/home/davidb/DESTv2_data_paper/15.demographic_inference/data/sfss/"
    # These SFSs respectively encode variation in the mainland region, the Americas
    # region (with the mainland region collapsed into one population), and Europe,
    # so Fsts less than 0.15 for each of them support the parsimonious conclusions
    # made earlier with the model selection performed with demographic inference
    # with moments followed by collapsed log-likelihood-based testing.
    sfs_names = ["k8_mainland_p2.npy",
                 "k8_Americas_p2.npy",
                 "k8_Europe.npy"]
    sfs_files = [data_dir + name for name in sfs_names]
    sfss = [moments.Spectrum(np.load(sfs_file)) for sfs_file in sfs_files]

    get_Fsts(sfss)

def get_Fsts(sfs_names, sfss):
    for sfs, sfs_name in zip(sfss, sfs_names):
        print(sfs_name, "\t", sfs.Fst())

if __name__ == "__main__":
    main()
