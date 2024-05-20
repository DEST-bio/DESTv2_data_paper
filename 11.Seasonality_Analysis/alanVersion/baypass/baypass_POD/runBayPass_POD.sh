#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH -A berglandlab
#SBATCH -c 16
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem 9G
#SBATCH -o /scratch/aob2x/DEST2_analysis/seasonality/logs/dest_pod.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DEST2_analysis/seasonality/logs/dest_pod.%A_%a.err # Standard error

### sbatch --array=102 ~/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/baypass/baypass_conrastPerm/contrast_permutation.sh
### sacct -j 58525234

# ijob -A berglandlab_standard -c16 -p standard --mem=9G

baypass="/home/aob2x/baypass_afton/sources/g_baypass"
module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
# SLURM_ARRAY_TASK_ID=51

### generate simulated data
  Rscript /scratch/aob2x/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/baypass/baypass_POD/runBayPass_POD.sh --args ${SLURM_ARRAY_TASK_ID}

### translate jobId into slice
  subPool=$( echo $((${SLURM_ARRAY_TASK_ID} % 50)) )
  echo $subPool

### RUN POD
  cd /scratch/aob2x/dest2_baypass/pods_v2
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v2/anapod_${SLURM_ARRAY_TASK_ID}_1 \
  -nthreads 16 \
  -omegafile /scratch/aob2x/dest2_baypass/pods_v2/omega.subpod_${SLURM_ARRAY_TASK_ID} \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}

  cd /scratch/aob2x/dest2_baypass/pods_v2
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v2/anapod_${SLURM_ARRAY_TASK_ID}_2 \
  -nthreads 16 \
  -omegafile /scratch/aob2x/dest2_baypass/pods_v2/omega.subpod_${SLURM_ARRAY_TASK_ID} \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}

  cd /scratch/aob2x/dest2_baypass/pods_v2
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v2/anapod_${SLURM_ARRAY_TASK_ID}_3 \
  -nthreads 16 \
  -omegafile /scratch/aob2x/dest2_baypass/pods_v2/omega.subpod_${SLURM_ARRAY_TASK_ID} \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}

  cd /scratch/aob2x/dest2_baypass/pods_v2
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v2/anapod_${SLURM_ARRAY_TASK_ID}_4 \
  -nthreads 16 \
  -omegafile /scratch/aob2x/dest2_baypass/pods_v2/omega.subpod_${SLURM_ARRAY_TASK_ID} \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}

  cd /scratch/aob2x/dest2_baypass/pods_v2
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v2/anapod_${SLURM_ARRAY_TASK_ID}_5 \
  -nthreads 16 \
  -omegafile /scratch/aob2x/dest2_baypass/pods_v2/omega.subpod_${SLURM_ARRAY_TASK_ID} \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}
