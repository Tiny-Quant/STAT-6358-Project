#!/bin/bash

#SBATCH -A mmcgee_stat6358_1252_001_0001
#SBATCH --job-name=sample_run_scratch
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH -t 0-5:0:0
#SBATCH -o job_%A_%a.out          # %A = array job ID, %a = array task ID
#SBATCH -e job_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ataychameekiatchai@smu.edu
#SBATCH --array=0-99%10

echo "Running task ID: $SLURM_ARRAY_TASK_ID"

cd /lustre/work/client/users/ataychameekiatchai/STAT-6358-Project/data || exit 1
pwd

module use /users/ataychameekiatchai/STAT-6358-Project/.devcontainer/
module load devcontainer

singularity exec \
  --bind /users/ataychameekiatchai/STAT-6358-Project \
  --pwd /users/ataychameekiatchai/STAT-6358-Project/data/gen_samples \
  /users/ataychameekiatchai/STAT-6358-Project/.devcontainer/devcontainer_latest_amd64_2025_04_17_13_08_09.sif \
  bash generate_sample.sh "$SLURM_ARRAY_TASK_ID"

pwd

# END OF SCRIPT
