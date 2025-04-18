#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -A berglandlab
#SBATCH -c 16
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem 9G
#SBATCH -o /scratch/aob2x/DEST2_analysis/seasonality/logs/dest_pod.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DEST2_analysis/seasonality/logs/dest_pod.%A_%a.err # Standard error

### sbatch --array=1-500 ~/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/baypass/baypass_POD/runBayPass_POD.v3.sh
### sacct -j 5676089
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/dest_pod.5676089_409.err
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/dest_pod.5676089_409.out

# ijob -A berglandlab_standard -c16 -p standard --mem=9G

baypass="/home/aob2x/baypass_afton/sources/g_baypass"
module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
# SLURM_ARRAY_TASK_ID=102

### generate simulated data
  Rscript --vanilla ~/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/baypass/baypass_POD/generatePOD_input.R ${SLURM_ARRAY_TASK_ID}

### translate jobId into slice
  subPool=$( echo $((${SLURM_ARRAY_TASK_ID} % 50)) )
  echo $subPool

### RUN POD
  cd /scratch/aob2x/dest2_baypass/pods_v3
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v3/anapod_${SLURM_ARRAY_TASK_ID}_1 \
  -poolsizefile   /standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_"${subPool}".poolsize \
  -nthreads 16 \
  -omegafile /standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_${subPool}_1_mat_omega.out \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}

  cd /scratch/aob2x/dest2_baypass/pods_v3
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v3/anapod_${SLURM_ARRAY_TASK_ID}_2 \
  -poolsizefile   /standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_"${subPool}".poolsize \
  -nthreads 16 \
  -omegafile /standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_${subPool}_2_mat_omega.out \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}

  cd /scratch/aob2x/dest2_baypass/pods_v3
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v3/anapod_${SLURM_ARRAY_TASK_ID}_3 \
  -poolsizefile   /standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_"${subPool}".poolsize \
  -nthreads 16 \
  -omegafile /standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_${subPool}_3_mat_omega.out \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}

  cd /scratch/aob2x/dest2_baypass/pods_v3
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v3/anapod_${SLURM_ARRAY_TASK_ID}_4 \
  -poolsizefile   /standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_"${subPool}".poolsize \
  -nthreads 16 \
  -omegafile /standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_${subPool}_4_mat_omega.out \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}

  cd /scratch/aob2x/dest2_baypass/pods_v3
  $baypass \
  -gfile /scratch/aob2x/dest2_baypass/pods_v2/G.subpod_${SLURM_ARRAY_TASK_ID} \
  -outprefix /scratch/aob2x/dest2_baypass/pods_v3/anapod_${SLURM_ARRAY_TASK_ID}_5 \
  -poolsizefile   /standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_"${subPool}".poolsize \
  -nthreads 16 \
  -omegafile /standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_${subPool}_5_mat_omega.out \
  -contrastfile /standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt \
  -seed ${SLURM_ARRAY_TASK_ID}
