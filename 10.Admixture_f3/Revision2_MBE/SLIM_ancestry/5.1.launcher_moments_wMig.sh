#!/bin/sh
#
#SBATCH -J moMIG
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH -o ./slurmOut/admix.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-999

module load python3.12-anaconda/2024.06-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda activate moments_jcbn
moment_script=/gpfs2/scratch/jcnunez/DEST2.0_analysis/REVISION3_ADMIX/5.SLIM_moments_admixWMIG.py

for k in {2..8}; do
echo $k

  python $moment_script \
	${SLURM_ARRAY_TASK_ID} \
	$k
	
done

conda deactivate

### Print the time
  echo "ended at"  `date`
