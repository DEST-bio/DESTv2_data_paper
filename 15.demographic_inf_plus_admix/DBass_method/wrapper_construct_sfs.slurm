#!/bin/bash

#SBATCH --mem 10G
#SBATCH --time 6:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab_standard
#SBATCH -c 1
#SBATCH -e slurm/wcs.err
#SBATCH --output slurm/wcs.out

# Command line argument is `popinfo_name_no_ext`
python construct_sfs.py $1
