#!/bin/bash
#SBATCH --job-name=bci-real-hybrid2
#SBATCH --output=logs/bci-real-hybrid2_%A_%a.out
#SBATCH --error=logs/bci-real-hybrid2_%A_%a.err
#SBATCH --time=24:00:00            # Time limit (HH:MM:SS)
#SBATCH --ntasks=1                 # Number of tasks (1 task for serial R scripts)
#SBATCH --cpus-per-task=16          # Number of CPU cores per task
#SBATCH --partition=gpu_4,gpu_8,gpu_4_a100,gpu_4_h100       # Partition to use (adjust to your system)
#SBATCH --mail-type=ALL            # Notifications for job start, end, fail
#SBATCH --mail-user=yannek.kaeber@biom.uni-freiburg.de # Your email for notifications
#SBATCH --mem=100G
#SBATCH --gres=gpu:1

################################################################################
# USAGE EXAMPLE:
#   sbatch run_torch_script.sh /path/to/your_script.R "1-10"
#
# The second argument (optional) specifies job indexes for the array.
################################################################################

# 1) Check that at least one argument (R script path) is provided
if [ "$#" -lt 1 ]; then
  echo "Usage: sbatch $0 <path_to_your_R_script> [array_indexes]"
  exit 1
fi

RSCRIPT_PATH="$1"

# 2) Check if a second argument is provided (defines Slurm array indexes)
if [ -n "$2" ]; then
  JOB_ARRAY_INDEXES="$2"
else
  JOB_ARRAY_INDEXES="1-32"  # Default array range
fi

# 3) Set Slurm array dynamically
#SBATCH --array=${JOB_ARRAY_INDEXES}

# 4) Load Singularity/Apptainer (adjust module as needed)
#module load system/singularity/3.11.3

# 5) (Optional) Set environment variables or simpler locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export R_PROFILE_USER=/dev/null

# 6) Run the container with GPU support, passing the job array index to the R script
singularity exec --nv ~/finn-r-torch.sif Rscript "$RSCRIPT_PATH" "$SLURM_ARRAY_TASK_ID"

# (Optional) For interactive usage (commented out):
# singularity exec --nv ~/finn-r-torch.sif R
