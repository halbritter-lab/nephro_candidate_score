#!/bin/bash
#
# Identify project and set log path
#SBATCH --job-name=extract_variant_features
#SBATCH --output=slurm_out/slurm_%j.out

#
# Set a required running time for the job.
#SBATCH --time=7-00
#
# Reserve resouces in partition
#SBATCH --partition medium
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G

echo "date: $(date)"
echo "host: $(hostname)"
echo "conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"


python extract_variant_features.py 2>&1 | tee "slurm_out/verbose_out_$SLURM_JOB_ID.txt"

echo "end"


