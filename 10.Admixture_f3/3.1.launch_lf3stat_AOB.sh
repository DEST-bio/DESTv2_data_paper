#!/usr/bin/env bash
#
#SBATCH -J glmOmn # A single job name for the array
#SBATCH -c 7
#SBATCH -N 1 # on one node
#SBATCH -t 00:60:00 #<= this may depend on your resources
#SBATCH --mem 63G #<= this may depend on your resources
#SBATCH -p standard
#SBATCH -A berglandlab_standard
#SBATCH -o /scratch/aob2x/DEST2_analysis/seasonality/logs/f3.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DEST2_analysis/seasonality/logs/f3.%A_%a.err # Standard error

### sbatch --array=1-4912 /home/aob2x/DESTv2_data_paper/10.Admixture_f3/3.1.launch_lf3stat_AOB.sh
### sacct -j 59068285
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.58265536*.err
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/f3.59068285_21.err


module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

Rscript \
--vanilla \
~/DESTv2_data_paper/10.Admixture_f3/3.0.calculate.f3stats_AOB.R \
${SLURM_ARRAY_TASK_ID}
