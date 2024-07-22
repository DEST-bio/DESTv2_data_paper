import moments
import numpy as np

def main():
    data_dir = "/home/davidb/DESTv2_data_paper/15.demographic_inference/data/sfss/"
    sfs_names = ["k8_mainland_p2.npy",
                 "k8_Americas_p2.npy",
                 "k8_Europe.npy"]
    sfs_files = [data_dir + name for name in sfs_names]
    sfss = [moments.Spectrum(np.load(sfs_file)) for sfs_file in sfs_files]

    get_Fsts(sfss)

def get_Fsts(sfss):
    for sfs in sfss:
        print(sfs.Fst())

if __name__ == "__main__":
    main()
