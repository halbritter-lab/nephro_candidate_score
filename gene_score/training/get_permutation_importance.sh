#!/bin/bash
#
# Identify project and set log path
#SBATCH --job-name=permutation_importance
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


# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --ID)
            MY_ID="$2"
            shift # past argument
            shift # past value
            ;;
        *)  
            # Unknown option
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if the required argument is provided
if [ -z "$MY_ID" ]; then
    echo "Please provide the --myID argument."
    exit 1
fi

# # Your existing script logic goes here, and you can use $MY_ID as needed
# echo "Executing script with MY_ID=$MY_ID"


file_path="/data/gpfs-1/users/rankn_c/work/halbritter/nephro_candidate_score/gene_score/training/get_permutation_importance.py"

# python "$file_path" 2>&1 | tee slurm_out/verbose_out_$SLURM_JOB_ID.txt

python "$file_path" --ID "$MY_ID" 2>&1 | tee slurm_out/verbose_out_$SLURM_JOB_ID.txt

echo "end"

# Example: sbatch command
# sbatch your_sbatch_script.sh --myID $MY_ID


