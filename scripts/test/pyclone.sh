#!/bin/bash
#SBATCH --job-name=PyClone
#SBATCH --partition=amd
#SBATCH --time=72:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.out
#SBATCH --error=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

## "==========================================================================="
## Credited from PyClone (https://github.com/Roth-Lab/pyclone)
## Authors: Zhong Guorui
## Date: 2025-09-05
## "==========================================================================="

## Activate the conda environment
conda activate pyclone

export MPLBACKEND=Agg
export PYCLONE_DATA_DIR="/lustre1/g/path_my/pipeline/somatic_variants_calling/data/DFSP/PyClone"

## Run PyClone for each case
case_ids=$(find "${PYCLONE_DATA_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)

## Run PyClone for each case
for case_id in $case_ids; do

    case_dir="${PYCLONE_DATA_DIR}/${case_id}"

    echo "$(date +"%F") $(date +"%T") - ${case_id} - Running PyClone ... "    
    # Change to the case directory
    cd "${case_dir}" || {
        echo "Error: Cannot change to directory ${case_dir}"
        continue
    }

    tsv_files=$(ls *.tsv)
    # echo "Files: $tsv_files"

    PyClone run_analysis_pipeline --in_files ${tsv_files} --working_dir pyclone_analysis >& pyclone.log

done
