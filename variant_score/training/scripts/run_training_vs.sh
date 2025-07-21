#!/bin/bash

# SLURM settings
#SBATCH --job-name=ML_training_vs
#SBATCH --output=slurm_%j.out  # Temporary; tee will handle the real output logging
#SBATCH --time=7-00
#SBATCH --partition medium
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G

echo "date: $(date)"
echo "host: $(hostname)"
echo "conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Assume CONFIG_YAML is set in the environment, e.g.:
# export CONFIG_YAML=/path/to/config.yaml
ML_PROJECTS_DIR=$(yq '.nephrology.ML_projectsdir' "$CONFIG_YAML")

SCRIPT_PATH="$ML_PROJECTS_DIR/variant_score/training/scripts/run_training_vs.py"
LOG_DIR="$ML_PROJECTS_DIR/variant_score/training/slurm_out"

mkdir -p "$LOG_DIR"  # Ensure log directory exists

python "$SCRIPT_PATH" 2>&1 | tee "$LOG_DIR/verbose_out_$SLURM_JOB_ID.txt"

echo "end"

