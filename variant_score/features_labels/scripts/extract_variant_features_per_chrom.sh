#!/bin/bash
#
# Identify project and set log path
#SBATCH --job-name=extract_variant_features_per_chrom
#SBATCH --output=slurm_out/slurm_%j.out

#
# Set a required running time for the job.
#SBATCH --time=4-00
#
# Reserve resouces in partition
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
file_path="${project_path}/nephro_candidate_score/variant_score/features_labels/scripts/extract_variant_features_per_chrom.py"


# parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --CHROM)
            chrom="$2"
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

# check if the required CHROM argument is provided
if [ -z "$chrom" ]; then
    echo "Please provide the --CHROM argument."
    exit 1
fi





# run python script
python "$file_path" --CHROM "$chrom" 2>&1 | tee "${project_path}/nephro_candidate_score/variant_score/features_labels/slurm_out/verbose_out_chr${chrom}_${SLURM_JOB_ID}.txt"


echo "end"



