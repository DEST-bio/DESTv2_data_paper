#!/usr/bin/env bash
#
#SBATCH -J glmOmn # A single job name for the array
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 03:30:00 #<= this may depend on your resources
#SBATCH --mem 15G #<= this may depend on your resources
#SBATCH -p standard
#SBATCH -A berglandlab
#SBATCH -o /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.%A_%a.err # Standard error

### sbatch --array=1-2173 /home/aob2x/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/launch.seasonality.Europe_29_sept_2023.sh
### sbatch --array=19,22,23,24,25,26,36,39,40,51,69,71,74,78,85,95,118,126,145,199,203,218,224,243,244,269,270,299,328,343,345,396,397,398,399,400,722,736,740,748,798,827,845,853,854,862,887,888,889,890,891,903,906,909,1005,1033,1052,1063,1066,1101,1105,1144,1146,1186,1190,1191,1192,1194,1206,1211,1214,1221,1230,1247,1811,1812,1813,1841,1850,1897,1900,1910,1935,1942,1966,1997,2022,2071,2086,2124,2127,2133,2134,2137,2138,2145
### sacct -j 53671778 --format="JobID%30,JobName,ExitCode,State" | grep "TIMEOUT" | sed -e 's/\s\+/,/g' | cut -f2 -d',' | cut -f2 -d'_' | tr '\n' ','
### seff 50160643_1
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.50142438_1000.out
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.53669286_1.err



module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1

cd /home/aob2x/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion


Rscript \
--vanilla \
seasonality_v2.R \
${SLURM_ARRAY_TASK_ID} \
"NoCore20_seas_europe" \
100
