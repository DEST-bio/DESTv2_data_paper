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

### i=1;k=str(2)
i=int(sys.argv[1])
k=str(sys.argv[2])
print("now processing", i, k)
###

Ancestral_id1="A_parent1"
Ancestral_id2="B_parent2"
Derived_id="C_derived"
root_dadi="dadi_sim_objects"
namsam = ["probs","0", "1", str(k),"rep",str(i),"delim"]
paths = [root_dadi,".".join(namsam)]
"/".join(paths)
fs_file = "/".join(paths)
###
pair_base = ["pop",k,"rep",str(i)]
Pair_name = ".".join(pair_base)
###
### collect the metadata
root_meta="./L_meta_sim_objects"
namsam2 = ["probs","0", "1", str(k),"rep",str(i),"meta"]
paths2 = [root_meta,".".join(namsam2)]
"/".join(paths2)
meta_file = "/".join(paths2)
inmet = pd.read_csv(meta_file,header=0,delimiter="\t")
Ancestral_n1=int(inmet.an1_ne[0])
Ancestral_n2=int(inmet.an2_ne[0])
Derived_n=int(inmet.Der_ne[0])
L = inmet.L[0]

print("fs file is", fs_file )
print("Pair_name is", Pair_name )
print("names are", Ancestral_id1,Ancestral_id2,Derived_id )

#constants
mu = 1.5e-6 #from SLiM
g = 1 
iterations=1
##### iterations=1

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

####
#opening output file to give column names
PMmod=open('./moments_slim_output/%s_output.admix.txt' % Pair_name,'w')
PMmod.write(
            str("Pair_name")+'\t'+ #print pair name
            str("L")+'\t'+ #double checking L is working as I want it to
            str("admix_prop")+'\t'+ #nu1
            str("fs_sanitycheck")+'\t'+
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()


print('defining functions')
####
# For modeling DEST data, pop_ids=[Afr, EU, NA], with Afr and EU interchangeable,
# so that NA is described as the result of an African-European admixture event.
def get_mig_mat(n_pops: int, m: float) -> np.ndarray:
    """
    :param n_pops int: Number of populations
    :param m float: Migration rate, in coalescent units of `1 / 2N_anc`
    :return: Square `n_pops` x `n_pops` matrix representing symmetric 
        migration, equals `m` everywhere except for zeros along the diagonal
    :rtype: np.ndarray
    """
    return m * np.ones([n_pops] * 2) - np.diag([m] * n_pops)

def admixture(params, ns, pop_ids=None):
    nu1, nu2, nu_admix, T_split, T_admix, mAE, mEA, admix_prop = params
    #mig_mat2 = get_mig_mat(2, m2)
    #mig_mat3 = get_mig_mat(3, m3)
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + 2 * ns[2])
    fs = moments.Spectrum(sts)
    fs = fs.split(0, ns[0] + ns[2], ns[1] + ns[2])
    mig_mat2 = [
    [0,   mAE],
    [mEA, 0  ]]
    fs.integrate([nu1, nu2], T_split, m=mig_mat2)
    fs = fs.admix(0, 1, ns[2], admix_prop)
    mig_mat3 = [
    [0,   mAE, 0],
    [mEA, 0  , 0],
    [0,   0  , 0]]
    fs.integrate([nu1, nu2, nu_admix], T_admix, m=mig_mat3)
    fs.pop_ids = pop_ids
    return fs.fold()
###   idx0, idx1, num_lineages, proportion, new_id=None
    
func_moments=admixture

params_guess = [
1.0, #nu1
1.0, #nu2
1.0, #nu_admix
0.5, #T_split
0.1, #T_admix
0.01, #m2
0.01, #m3
0.5 # admix_prop
]

lower_bounds = [
0.01,   # nu1
0.01,   # nu2
0.01,   # nu_admix
0.001,  # T_split
0.001,  # T_admix
0.001,  #m2
0.001,  #m3
0.01 # admix_prop
]

upper_bounds = [
10.0,   # nu1
10.0,   # nu2
10.0,   # nu3
10.0,   # T_split
5.0,    # T_admix
0.1,  #m2
0.1,  #m3
0.99 # admix_prop
]
#upper_bound = [100, 100, 100, 100, 100, 1]
#lower_bound = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 0]


print('optimization loop')
for i in range(int(iterations)): #iterations is imported from sys. argument #1
	print("starting optimization "+str(i))
#This number is 5. i.e., count parameters. 
	params = len(["nu1", "nu2", "nu_admix", "T_split", "T_admix", "mAE", "mEA", "admix_prop"]) #for use in AIC calculation
	#Start the run by picking random parameters from a uniform distribution.
	#popt=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(params)]
	    
	#This is the optimization step for moments.
	#popt is the prior.
	#fs folded is a tranform SFS by folding it. The original SFS loaded in sys.arg #1 is quasi-folded. i.e. polarized to reference genome.
	#Folding it is a must, because reference is not 100% ancestral
	popt=moments.Inference.optimize_log(params_guess, fs_folded, 
	func_moments,
	lower_bound=lower_bounds, 
	upper_bound=upper_bounds,
	verbose=True, 
	maxiter=100)
	
	model = func_moments(popt, ns)
	
	#Calculate log likelihood of the model fit
	ll_model=moments.Inference.ll_multinom(model, fs_folded)
	#Now calculate AIC of model fit
	aic = 2*params - 2*ll_model
	print('Maximum log composite likelihood: {0}'.format(ll_model))
	
	#Now estimate theta from model fit
	theta = moments.Inference.optimal_sfs_scaling(model, fs_folded)
	#reducing complexity of calculations to follow by adding variables in lieu of expressions/esoteric df calls
	Nref= theta/(4*mu*L)
	nu1=popt[0]
	nu2=popt[1]
	nu3=popt[2] 
	T_split=popt[3] 
	T_admix=popt[4]  
	admix_prop=popt[7] ## in this model this is position 7!
	
	#Open the output file
	PMmod=open('./moments_slim_output/%s_output.admix.txt' % Pair_name,'a')
	    #Dumping output ot outfile
	PMmod.write(
	        str(Pair_name)+'\t'+ #print pair name
	        str(fs_file)+'\t'+ #double checking fs is the right one
	        str(L)+'\t'+ #double checking L is working as desired
	        str(admix_prop)+'\t'+ #nu1
	        str(S)+'\t'+ #sanity check... should give number of segregating sites in SFS
	        str(ll_model)+'\t'+
	        str(aic)+'\n')
	PMmod.close()
	
	print("Moments finished running")
