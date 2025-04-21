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

for i in range(114): #iterations is imported from sys. argument #1
	print("starting pickling "+str(i))
	indat = pd.read_csv("pairs_guide_file.txt",header=0,delimiter="\t")
	Ancestral_id1=indat.Parent1[i]
	Ancestral_id2=indat.Parent2[i]
	Derived_id=indat.Derived[i]
	root_dadi="dadi_objects"
	namsam = ["probs",Ancestral_id1, Ancestral_id2, Derived_id,"delim"]
	paths = [root_dadi,".".join(namsam)]
	"/".join(paths)
	fs_file = "/".join(paths)
	Pair_name = Derived_id
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
	dd = dadi.Misc.make_data_dict(fs_file) #reads in genomalicious SNP file
	data = pd.read_csv(fs_file, sep="\t", nrows=1)
	pop_id=[Ancestral_id1,Ancestral_id2,Derived_id]
	pools=[Ancestral_n1, Ancestral_n2,Derived_n ] 
	fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=pools, polarized=False) #takes data dict and folds
	with open('sfs_objs/%s.sfs.plk' % Pair_name, 'wb') as f:
		pickle.dump(fs_folded, f, pickle.HIGHEST_PROTOCOL)

#### Note to self: Here I am taking a pickling strategy due to an architectural issue with the VACC... 