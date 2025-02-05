#!/bin/bash

# Identify project and set log path
#SBATCH --job-name=ML_training
#SBATCH --output=gene_score/training/slurm_out/slurm_%j.out

# Set a required running time for the job.
#SBATCH --time=7-00

# Reserve resouces in partition
#SBATCH --partition medium
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G

echo "date: $(date)"
echo "host: $(hostname)"
echo "conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# # Run the command and capture its output
# project_path=$(yq -r .nephrology.ML_projectsdir "$CONFIG_FILE")

# # Now $file_path contains the result
# # echo "The file path is: $project_path"

# # append project path
# gs_training_path="${project_path}nephro_candidate_score/gene_score/training"


# # gs_training_path="$project_path/nephro_candidate_score/gene_score/training"

# python "$gs_training_path/run_training.py" 2>&1 | tee $gs_training_path/slurm_out/verbose_out_$SLURM_JOB_ID.txt


file_path="/data/gpfs-1/users/rankn_c/work/halbritter/nephro_candidate_score/gene_score/training/scripts/run_training.py"

python "$file_path" 2>&1 | tee gene_score/training/slurm_out/verbose_out_$SLURM_JOB_ID.txt


echo "end"
