#!/usr/bin/env bash
#
#SBATCH -J glmOmn # A single job name for the array
#SBATCH -c 20
#SBATCH -N 1 # on one node
#SBATCH -t 00:30:00 #<= this may depend on your resources
#SBATCH --mem 150G #<= this may depend on your resources
#SBATCH -p standard
#SBATCH -A berglandlab_standard
#SBATCH -o /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.%A_%a.err # Standard error

### sbatch --array=1-101 /home/aob2x/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/windowEnrichmentTests/launch.windowAnalysis_BAYPASS_CONTRAST_perm.sh
### sacct -j 5675159 --format="JobID%30,JobName,ExitCode,State" | grep "TIMEOUT" | sed -e 's/\s\+/,/g' | cut -f2 -d',' | cut -f2 -d'_' | tr '\n' ','
### seff 5675159_1
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.5675159_10.err
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.53669286_1.err

module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

Rscript \
--vanilla \
~/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/windowEnrichmentTests/windowAnalysis_BAYPASS_CONTRAST_perm.R \
${SLURM_ARRAY_TASK_ID}
