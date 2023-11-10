#!/bin/bash
#
# Identify project and set log path
#SBATCH --job-name=ML_training
#SBATCH --output=slurm_out/slurm_%j.out

#
# Set a required running time for the job.
#SBATCH --time=1-00
#
# Reserve resouces in partition
#SBATCH --partition medium
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G

echo "date: $(date)"
echo "host: $(hostname)"
echo "conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

file_path="/data/gpfs-1/users/rankn_c/work/halbritter/nephro_candidate_score/gene_score/training/run_training.py"

python "$file_path" 2>&1 | tee slurm_out/verbose_out_$SLURM_JOB_ID.txt

echo "end"
