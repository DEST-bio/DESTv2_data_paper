#!/bin/bash
#SBATCH -J dest_subpool # A single job name for the array
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 9G
#SBATCH -t 20:30:00
#SBATCH -p standard
#SBATCH -A berglandlab_standard
#SBATCH -o /scratch/aob2x/logs/dest_baypass.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/dest_baypass.%A_%a.err # Standard error

### sbatch --array=1-2500 ~/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/baypass/runBaypass/dest_subpoolsarray_perm.sh
### sacct -j 58264320
### cat /scratch/aob2x/logs/dest_baypass.58264320_1.out

# ijob -A berglandlab -c16 -p standard --mem=10G
# SLURM_ARRAY_TASK_ID=1

# head /standard/vol186/bergland-lab/alan/dest_baypass/subpoolcontrol.txt

a=1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /standard/vol186/bergland-lab/alan/dest_baypass/subpoolcontrol_perm.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo $opt2
echo $opt1
echo $opt3

module load gcc/11.4.0

baypass="/home/aob2x/baypass_public/sources/g_baypass"

cd /standard/vol186/bergland-lab/alan/dest_baypass/dest_permoutput

if [ ! -f /standard/vol186/bergland-lab/Gio/dest2_season_contrast_"${opt3}"_long.txt ]; then
  cat /standard/vol186/bergland-lab/Gio/dest2_season_contrast_"${opt3}".txt | tr '\n' ' ' | sed 's/ $/\n/g' > /standard/vol186/bergland-lab/Gio/dest2_season_contrast_"${opt3}"_long.txt
fi


$baypass -gfile /standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_"${opt1}".genobaypass \
-poolsizefile /standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_"${opt1}".poolsize \
-omegafile  /standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_"${opt1}"_"${opt2}"_mat_omega.out \
-outprefix "destsubpool_${opt1}_${opt2}_${opt3}" \
-nthreads 16  \
-contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast_"${opt3}"_long.txt \
-seed $opt2
