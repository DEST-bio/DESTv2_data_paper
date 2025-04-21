#!/bin/sh
#
#SBATCH -J moments_admix
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH -o ./slurmOut/admix.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-114

module load python3.12-anaconda/2024.06-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda activate moments_jcbn


moment_script=/gpfs2/scratch/jcnunez/DEST2.0_analysis/MOMENTS_REVISION/3.run_moments_admix.py

guide=/gpfs2/scratch/jcnunez/DEST2.0_analysis/MOMENTS_REVISION/pairs_guide_file.txt

Ancestral_id1=$( cat $guide | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F "\t" '{ print $1 }' )
Ancestral_id2=$( cat $guide | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F "\t" '{ print $2 }' )
Derived_id=$( cat $guide | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F "\t" '{ print $3 }' )

echo $Ancestral_id1 $Ancestral_id2 $Derived_id

guide_meta=/gpfs2/scratch/jcnunez/DEST2.0_analysis/MOMENTS_REVISION/L_meta_objects/probs.${Ancestral_id1}.${Ancestral_id2}.${Derived_id}.meta

Ancestral_n1=$( cat $guide_meta | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F "\t" '{ print $1 }' )
Ancestral_n2=$( cat $guide_meta | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F "\t" '{ print $2 }' )
Derived_n=$( cat $guide_meta | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F "\t" '{ print $3 }' )
L=$( cat $guide_meta | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F "\t" '{ print $4 }' )

echo $Ancestral_n1 $Ancestral_n2 $Derived_n $L

Pair_name=${Derived_id}
SFS=/gpfs2/scratch/jcnunez/DEST2.0_analysis/MOMENTS_REVISION/dadi_objects/probs.${Ancestral_id1}.${Ancestral_id2}.${Derived_id}.delim

####
### Finally a sanity check
  echo "Now Processing" $Pair
  echo "Now Loading" $SFS "=> where" $Pair "SFS is located"
  echo $Pair "has an L parameter of" $L "bp"

 mkdir "est_output"

 cd ./est_output

  python $moment_script \
	${SFS} \
	${Pair_name} \
    ${Ancestral_id1} \
    ${Ancestral_id2} \
    ${Derived_id} \
    ${Ancestral_n1} \
    ${Ancestral_n2} \
    ${Derived_n} \
    10 \
    ${L}

conda deactivate

### Print the time
  echo "ended at"  `date`
