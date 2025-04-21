#module load python3.12-anaconda/2024.06-1;conda activate moments_jcbn; python

# import packages that'll be used
import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
import dadi
from dadi import Misc
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from datetime import datetime

###
i=int(sys.argv[1])
print("now processing", i)
###

indat = pd.read_csv("pairs_guide_file.txt",header=0,delimiter="\t")
Ancestral_id1=indat.Parent1[i]
Ancestral_id2=indat.Parent2[i]
Derived_id=indat.Derived[i]
root_dadi="dadi_objects"
namsam = ["probs",Ancestral_id1, Ancestral_id2, Derived_id,"delim"]
paths = [root_dadi,".".join(namsam)]
"/".join(paths)
fs_file = "/".join(paths)
###
Pair_name = Derived_id
###
### collect the metadata
root_meta="./L_meta_objects"
namsam2 = ["probs",Ancestral_id1, Ancestral_id2, Derived_id,"meta"]
paths2 = [root_meta,".".join(namsam2)]
"/".join(paths2)
meta_file = "/".join(paths2)
inmet = pd.read_csv(meta_file,header=0,delimiter="\t")
Ancestral_n1=inmet.an1_ne[0]
Ancestral_n2=inmet.an2_ne[0]
Derived_n=inmet.Der_ne[0]
L = inmet.L[0]

print("fs file is", fs_file )
print("Pair_name is", Pair_name )
print("names are", Ancestral_id1,Ancestral_id2,Derived_id )

#constants
mu = 2.8e-9 #from Keightley et al. 2014
g = 0.06666667 #equals 15 gen/year --- See Pool 2015
iterations=10
#####

####
print('loading data')
dd = dadi.Misc.make_data_dict(fs_file) #reads in genomalicious SNP file
data = pd.read_csv(fs_file, sep="\t", nrows=1)

#####
#setting pop id's and projections from if/else
pop_id=["A_parent1",  "B_parent2",  "C_derived"]
pools=[Ancestral_n1, Ancestral_n2,Derived_n ] 

print('folding sfs')
fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=pools, polarized=False) #takes data dict and folds
ns = fs_folded.sample_sizes #gets sample sizes of dataset
S = fs_folded.S()

def admixture(params, ns):
	nu1, nu2, nu_admix, T_split, T_admix, p_admix = params
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + 2 * ns[2])
	fs = moments.Spectrum(sts)
	fs = moments.Manips.split_1D_to_2D(fs, ns[0] + ns[2], ns[1] + ns[2])
	fs.integrate([nu1, nu2], T_split, m=np.zeros([2, 2]))
	fs = fs.admix(0, 1, ns[2], p_admix)
	fs.integrate([nu1, nu2, nu_admix], T_admix, m=np.zeros([3, 3]))
	return fs.fold()
    
# Initial parameter guesses: [T_div, T_admix, f]
params_guess = [
    1.0,   # nu1
    1.0,   # nu2
    1.0,   # nu_admix
    0.5,   # T_split
    0.1,   # T_admix
    0.5    # admix_prop (even mix)
]

lower_bounds = [
    0.01,   # nu1
    0.01,   # nu2
    0.01,   # nu_admix
    0.001,  # T_split
    0.001,  # T_admix
    0.01    # admix_prop (avoid 0/1 to keep numerical stability)
]

upper_bounds = [
    10.0,   # nu1
    10.0,   # nu2
    10.0,   # nu_admix
    10.0,   # T_split
    5.0,    # T_admix
    0.99    # admix_prop
]

# Optimize
popt=moments.Inference.optimize_log(
params_guess, 
fs_folded, 
admixture,
lower_bound=lower_bounds, 
upper_bound=upper_bounds,
verbose=True, 
maxiter=100)
