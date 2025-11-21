#!/bin/bash
#SBATCH --job-name=pyclone-vi
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
## Credited from pyclone-vi (https://github.com/Roth-Lab/pyclone-vi/tree/master)
## Authors: Zhong Guorui
## Date: 2025-09-05
## "==========================================================================="

## Create pyclone-vi conda environment
# conda create -c conda-forge -n pyclone-vi --file /lustre1/g/path_my/pipeline/somatic_variants_calling/conf/pyclone-vi.txt --yes
conda activate pyclone-vi
# pip install git+https://github.com/Roth-Lab/pyclone-vi.git
# pyclone-vi --help
# pyclone-vi fit -i /lustre1/g/path_my/pipeline/somatic_variants_calling/data/examples/tracerx.tsv -o tracerx.h5 -c 40 -d beta-binomial -r 10
# pyclone-vi write-results-file -i tracerx.h5 -o tracerx.tsv

export PYCLONE_DATA_DIR="/lustre1/g/path_my/pipeline/somatic_variants_calling/data/DFSP/pyclone-vi"

## Run PyClone for each case
case_ids=$(find "${PYCLONE_DATA_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)

## Run PyClone for each case
for case_id in $case_ids; do

    input_file="${PYCLONE_DATA_DIR}/${case_id}"/${case_id}.tsv
    output_file="${PYCLONE_DATA_DIR}/${case_id}"/${case_id}.h5
    
    echo "$(date +"%F") $(date +"%T") - ${case_id} - Running PyClone ... "    

    pyclone-vi fit \
        -i "${input_file}" \
        -o "${output_file}" \
        -c 40 \
        -d beta-binomial \
        -r 10
    
    pyclone-vi write-results-file -i tracerx.h5 -o tracerx.tsv


done
