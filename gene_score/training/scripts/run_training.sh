#!/bin/bash

#SBATCH --job-name=ML_training
#SBATCH --output=gene_score/training/slurm_out/slurm_%j.out
#SBATCH --time=7-00
#SBATCH --partition medium
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G

echo "date: $(date)"
echo "host: $(hostname)"
echo "conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Check if CONFIG_FILE is set
if [ -z "$CONFIG_FILE" ]; then
    echo "Error: CONFIG_FILE environment variable is not set."
    exit 1
fi

# Extract project path from YAML config
project_path=$(yq -r .nephrology.ML_projectsdir "$CONFIG_FILE")

# Construct script path
file_path="${project_path}/nephro_candidate_score/gene_score/training/scripts/run_training.py"

# Run the training script
python "$file_path" 2>&1 | tee gene_score/training/slurm_out/verbose_out_$SLURM_JOB_ID.txt

echo "end"




