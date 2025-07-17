#!/bin/bash
#
#SBATCH --job-name=permutation_importance
#SBATCH --output=gene_score/training/slurm_out/slurm_%j.out
#SBATCH --time=1-00
#SBATCH --partition medium
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G

echo "date: $(date)"
echo "host: $(hostname)"
echo "conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Ensure CONFIG_FILE is set
if [ -z "$CONFIG_FILE" ]; then
    echo "Error: CONFIG_FILE environment variable is not set."
    exit 1
fi

# Extract project path from config
project_path=$(yq -r .nephrology.ML_projectsdir "$CONFIG_FILE")

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --ID)
            training_ID="$2"
            shift
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if ID argument was provided
if [ -z "$training_ID" ]; then
    echo "Please provide the --ID argument."
    exit 1
fi

# Construct path to Python script dynamically
file_path="${project_path}/nephro_candidate_score/gene_score/training/scripts/get_permutation_importance.py"

# Run script with dynamic path and logging
python "$file_path" --ID "$training_ID" 2>&1 | tee gene_score/training/slurm_out/verbose_out_$SLURM_JOB_ID.txt

echo "end"

