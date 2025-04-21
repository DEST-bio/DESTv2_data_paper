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
import pickle

###
i=sys.argv[1]
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
objsids = ["./sfs_objs/",Derived_id,".sfs.","plk"]
objsids_fil="".join(objsids)

with open(objsids_fil, 'rb') as f:
	fs_folded = pickle.load(f)

####
####
ns = fs_folded.sample_sizes #gets sample sizes of dataset
S = fs_folded.S()

####
#opening output file to give column names
PMmod=open('%s_output.admix.txt' % Pair_name,'w')
PMmod.write(
            str("Pair_name")+'\t'+ #print pair name
            str("L")+'\t'+ #double checking L is working as I want it to
            str("admix_prop")+'\t'+ #nu1
            str("fs_sanitycheck")+'\t'+
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()


####
# For modeling DEST data, pop_ids=[Afr, EU, NA], with Afr and EU interchangeable,
# so that NA is described as the result of an African-European admixture event.
def admixture(params, ns, pop_ids=None):
    nu1, nu2, nu3, T_split, T_admix, admix_prop = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + 2 * ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0] + ns[2], ns[1] + ns[2])
    fs.integrate([nu1, nu2], T_split)
    fs = fs.admix(0, 1, ns[2], admix_prop)
    fs.integrate([nu1, nu2, nu3], T_admix)
    fs.pop_ids = pop_ids
    return fs
###   
    
func_moments = admixture

upper_bound = [100, 100, 100, 100, 100, 1]
lower_bound = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 0]


for j in range(int(iterations)): #iterations is imported from sys. argument #1
	print("starting optimization "+str(i))
#This number is 5. i.e., count parameters. 
	params = len(["nu1", "nu2", "nu3", "T_split", "T_admix", "admix_prop"]) #for use in AIC calculation
	#Start the run by picking random parameters from a uniform distribution.
	popt=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(params)]
	    
	#This is the optimization step for moments.
	#popt is the prior.
	#fs folded is a tranform SFS by folding it. The original SFS loaded in sys.arg #1 is quasi-folded. i.e. polarized to reference genome.
	#Folding it is a must, because reference is not 100% ancestral
	popt=moments.Inference.optimize_log(popt, fs_folded, func_moments,lower_bound=lower_bound, upper_bound=upper_bound,verbose=False, maxiter=100,)
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
	admix_prop=popt[5] 
	
	#Open the output file
	PMmod=open('%s_output.admix.txt' % Pair_name,'a')
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