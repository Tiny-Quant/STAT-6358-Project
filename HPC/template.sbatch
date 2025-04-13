#!/bin/bash



#SBATCH -A mmcgee_stat6358_1252_001_0001



# Name for the job that will be visible in the job queue and accounting tools.

#SBATCH --job-name NewJob



# Name of the SLURM partition that this job should run on.

#SBATCH --cpus-per-task=16       # Adjusted CPU allocation
#SBATCH --mem=128G               # Increased memory allocation



# Time limit for the job in the format Days-H:M:S

# A job that reaches its time limit will be cancelled.

# Specify an accurate time limit for efficient scheduling so your job runs promptly.

#SBATCH -t 0-1:0:0



# The standard output and errors from commands will be written to these files.

# %j in the filename will be replace with the job number when it is submitted.

#SBATCH -o job_%j.out

#SBATCH -e job_%j.err



# Send an email when the job status changes, to the specfied address.

#SBATCH --mail-type ALL

#SBATCH --mail-user ataychameekiatchai@smu.edu



# Log common failure points.

cd /lustre/work/client/users/ataychameekiatchai/STAT-6358-Project/data
pwd



module use /users/ataychameekiatchai/STAT-6358-Project/.devcontainer/
module load devcontainer
singularity exec \
  --bind /users/ataychameekiatchai/STAT-6358-Project \
  --pwd /users/ataychameekiatchai/STAT-6358-Project/data \
  /users/ataychameekiatchai/STAT-6358-Project/.devcontainer/devcontainer_latest_amd64_2025_04_10_18_43_05.sif \
  bash pull_SRA_and_fastq.sh
pwd



# END OF SCRIPT
