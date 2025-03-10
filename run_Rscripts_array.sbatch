#!/bin/bash
#SBATCH --job-name=bci-real-genus-full
#SBATCH --output=bci-real-genus-full_%A_%a.out
#SBATCH --error=bci-real-genus-full_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu-single
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yannek.kaeber@biom.uni-freiburg.de
#SBATCH --mem=50G
#SBATCH --gres=gpu:1
#SBATCH --array=1-9

################################################################################
# USAGE EXAMPLE:
#   sbatch run_torch_script.sh /path/to/your_script.R
#
# This script will run that R script inside the Singularity container.
################################################################################

# 1) Check that an R script path was provided
if [ "$#" -ne 1 ]; then
echo "Usage: sbatch $0 <path_to_your_R_script>"
exit 1
fi

RSCRIPT_PATH="$1"

# 2) Load Singularity/Apptainer (adjust module as needed)
module load system/singularity/3.11.3

# 3) (Optional) Set environment variables or simpler locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export R_PROFILE_USER=/dev/null

# 4) Run the container with GPU support, passing the job array index to the R script
singularity exec --nv ~/finn-r-torch.sif Rscript "$RSCRIPT_PATH" "$SLURM_ARRAY_TASK_ID"

# (Optional) For interactive usage (commented out):
# singularity exec --nv ~/finn-r-torch.sif R
