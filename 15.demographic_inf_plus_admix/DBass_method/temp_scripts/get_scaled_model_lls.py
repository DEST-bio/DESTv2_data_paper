import numpy as np
import moments
import os
import pandas as pd
import moments_models as mm

def main():
    # Setup
    data_dir = "/project/berglandlab/david_DEST2/"
    sfs_dir = data_dir + "sfss/"
    collected_op_file = "output/collected_output.tsv"

    sfss = {f.rsplit(".")[0]: moments.Spectrum(np.load(sfs_dir + f)).fold() 
            for f in os.listdir(sfs_dir) if f.endswith(".npy")}
    df = pd.read_csv(collected_op_file,
                     names=["ll", "coll_pop_ll"] + \
                        [f'est{i}' for i in range(8)],
                        index_col=[0, 1])
    
    print(df)
    print()

    # Europe
    # lines [0-7)
    Europe_admixture0_sfs = mm.admixture(list(df.loc[('Europe', 'admixture_0')])[2:],
                                         sfss['k8_Europe'].sample_sizes)
    Europe_admixture1_sfs = mm.admixture(list(df.loc[('Europe', 'admixture_1')])[2:],
                                         sfss['k8_Europe'].sample_sizes)
    Europe_admixture2_sfs = mm.admixture(list(df.loc[('Europe', 'admixture_2')])[2:],
                                         sfss['k8_Europe'].sample_sizes)
    Europe_twosplits0_sfs = mm.twosplits(list(df.loc[('Europe', 'twosplits_0')])[2:],
                                         sfss['k8_Europe'].sample_sizes)
    Europe_twosplits1_sfs = mm.twosplits(list(df.loc[('Europe', 'twosplits_1')])[2:],
                                         sfss['k8_Europe'].sample_sizes)
    Europe_twosplits2_sfs = mm.twosplits(list(df.loc[('Europe', 'twosplits_2')])[2:],
                                         sfss['k8_Europe'].sample_sizes)
    Europe_split_sfs = mm.split(list(df.loc[('Europe', 'split')])[2:6],
                                sfss['k4_Europe'].sample_sizes)

    print("### Europe ###", end="\n\n")
    print("Admixture")
    coll_ll(Europe_admixture0_sfs, sfss['k8_Europe'].swapaxes(0, 2), 1)
    coll_ll(Europe_admixture1_sfs, sfss['k8_Europe'].swapaxes(1, 2), 2)
    coll_ll(Europe_admixture2_sfs, sfss['k8_Europe'], 1)
    print("Twosplits")
    coll_ll(Europe_twosplits0_sfs, sfss['k8_Europe'], 1)
    coll_ll(Europe_twosplits1_sfs, sfss['k8_Europe'].swapaxes(0, 1), 0)
    coll_ll(Europe_twosplits2_sfs, sfss['k8_Europe'].swapaxes(0, 2), 1) 
    print("Split")
    coll_ll(Europe_split_sfs, sfss['k4_Europe'], -1)

    # mainland
    # lines [7, 9)
    mainland_split_sfs = mm.split(list(df.loc[('mainland', 'split')])[2:6],
                                  sfss['k8_mainland_p2'].sample_sizes)
    mainland_two_epoch_sfs = mm.two_epoch(list(df.loc[('mainland', 'two_epoch')])[2:4],
                                          sfss['k8_mainland_p1'].sample_sizes)
    
    print("### mainland ###", end="\n\n")
    print("Split")
    coll_ll(mainland_split_sfs, sfss['k8_mainland_p2'], 1)
    print("Two epoch")
    coll_ll(mainland_two_epoch_sfs, sfss['k8_mainland_p1'], -1)

    # Americas
    # lines [9, 11)
    Americas_split_sfs = mm.split(list(df.loc[('Americas', 'split')])[2:6],
                                  sfss['k8_Americas_p2'].sample_sizes)
    Americas_two_epoch_sfs = mm.two_epoch(list(df.loc[('Americas', 'two_epoch')])[2:4],
                                          sfss['k8_Americas_p1'].sample_sizes)
    
    print("### Americas ###", end="\n\n")
    print("Two epoch")
    coll_ll(Americas_split_sfs, sfss['k8_Americas_p2'], 1)
    print("Split")
    coll_ll(Americas_two_epoch_sfs, sfss['k8_Americas_p1'], -1)

    # Transatlantinc
    # lines [11, 14)
    Transatlantic_admixture0_sfs = mm.admixture(list(df.loc[('Transatlantic', 'admixture_0')])[2:],
                                                sfss['k8_Transatlantic'].sample_sizes[1:])
    Transatlantic_admixture1_sfs = mm.admixture(list(df.loc[('Transatlantic', 'admixture_1')])[2:],
                                                sfss['k8_Transatlantic'].sample_sizes[1:])
    Transatlantic_admixture2_sfs = mm.admixture(list(df.loc[('Transatlantic', 'admixture_2')])[2:],
                                                sfss['k8_Transatlantic'].sample_sizes[1:])
    
    print("### Transatlantic ###", end="\n\n")
    print("Caribbean and Guinea admixture")
    coll_ll(Transatlantic_admixture0_sfs, sfss['k8_Transatlantic'].swapaxes(2, 3).\
        marginalize([0]), -1)
    print("Europe and Guinea admixture")
    coll_ll(Transatlantic_admixture1_sfs, sfss['k8_Transatlantic'].swapaxes(2, 3).\
        marginalize([1]), -1)
    print("Europe and Caribbean admixture")
    coll_ll(Transatlantic_admixture2_sfs, sfss['k8_Transatlantic'].swapaxes(2, 3).\
        marginalize([2]), -1)

def coll_ll(model_sfs: moments.Spectrum, 
            data: moments.Spectrum,
            coll_axis: int) -> None:
    
    ll_norm = moments.Inference.ll_multinom(model_sfs, data / data.S())
    #ll_norm = moments.Inference.ll_multinom(model_sfs, data)

    print(ll_norm)
    
    if coll_axis >= 0:
        coll_model_sfs = model_sfs.marginalize([coll_axis])
        coll_data = data.marginalize([coll_axis])

        ll_coll_norm = moments.Inference.ll_multinom(coll_model_sfs, coll_data / coll_data.S())
        #ll_coll_norm = moments.Inference.ll_multinom(coll_model_sfs, coll_data)
        
        print(ll_coll_norm)
    
    print()

if __name__ == "__main__":
    main()
