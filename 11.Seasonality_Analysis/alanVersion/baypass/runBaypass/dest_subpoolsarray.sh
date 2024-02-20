#!/bin/bash
#SBATCH -J dest_subpool # A single job name for the array
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 9G
#SBATCH -t 20:30:00
#SBATCH -p standard
#SBATCH -A berglandlab
#SBATCH -o /scratch/aob2x/logs/dest_baypass.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/dest_baypass.%A_%a.err # Standard error

### sbatch --array=1-250 /scratch/aob2x/GioMazzeoDESTWork/alan/dest_subpoolsarray.sh
### sacct -j 54474211
### cat /scratch/aob2x/logs/dest_baypass.54474211_1.out

# SLURM_ARRAY_TASK_ID=1
a=1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /standard/vol186/bergland-lab/alan/dest_baypass/subpoolcontrol.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo $opt2
echo $opt1

module load gcc/11.4.0

baypass="/home/aob2x/baypass_public/sources/g_baypass"

cd /standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass

$baypass -gfile /standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_"${opt1}".genobaypass \
-poolsizefile   /standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_"${opt1}".poolsize \
-outprefix "destsubpool_${opt1}_${opt2}" \
-nthreads 16  \
-contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
-seed $opt2


tar -czvf ~/baypass_example.tar.gz \
/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_21_2_summary_yij_pij.out \
/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_21_2_summary_contrast.out \
/standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt
