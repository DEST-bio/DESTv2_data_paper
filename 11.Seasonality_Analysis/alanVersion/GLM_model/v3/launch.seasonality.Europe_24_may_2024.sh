#!/usr/bin/env bash
#
#SBATCH -J glmOmn # A single job name for the array
#SBATCH -c 20
#SBATCH -N 1 # on one node
#SBATCH -t 03:30:00 #<= this may depend on your resources
#SBATCH --mem 100G #<= this may depend on your resources
#SBATCH -p standard
#SBATCH -A berglandlab
#SBATCH -o /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.%A_%a.err # Standard error

### sbatch --array=1-2173 /home/aob2x/DESTv2_data_paper/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/GLM_model/v3/launch.seasonality.Europe_24_may_2024.sh
###
### sacct -j 53885553 --format="JobID%30,JobName,ExitCode,State" | grep "TIMEOUT" | sed -e 's/\s\+/,/g' | cut -f2 -d',' | cut -f2 -d'_' | tr '\n' ','
### seff 50160643_1
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.50142438_1000.out
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.53669286_1.err



module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1

cd DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/GLM_model/v3/


Rscript \
--vanilla \
seasonality_v3.R \
${SLURM_ARRAY_TASK_ID} \
"NoCore20_seas_europe" \
100
