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
python features_labels/scripts/extract_variant_features_per_chrom.py --CHROM "$chrom" 2>&1 | tee features_labels/slurm_out/verbose_out_chr"$chrom"_$SLURM_JOB_ID.txt


echo "end"



