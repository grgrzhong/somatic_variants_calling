#!/bin/bash
#SBATCH --job-name=telomerehunter
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

source $(conda info --base)/etc/profile.d/conda.sh
## Create the environment for telomerehunter
## * Can only run in python 2.7 environment
# conda create -n telomerehunter python=2.7
conda activate telomerehunter
# pip install --no-cache-dir telomerehunter
# conda install r-reshape2 r-tidyverse
module load gnuparallel/20211222

export BAM_DIR=/lustre1/g/path_my/pipeline/somatic_variants_calling/data/DFSP/BAM
export TELOMERHUNTER_DIR=/lustre1/g/path_my/pipeline/somatic_variants_calling/data/DFSP/TelomereHunter
mkdir -p ${TELOMERHUNTER_DIR}

export BANDING_FILE="/lustre1/g/path_my/pipeline/somatic_variants_calling/data/Reference/hg38_banding.tsv"

## Function to run telomerehunter for a sample
run_telomerhunter() {
    
    local tumour_id="$1"

    patient_id=$(echo "$tumour_id" | cut -d'-' -f1,2)
    normal_id=${patient_id}-N
    

    if [ -d "${BAM_DIR}/${normal_id}" ]; then

        out_dir=${TELOMERHUNTER_DIR}/${tumour_id}
        mkdir -p "$out_dir"

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Running TelomereHunter ..."
    
        tumour_bam="${BAM_DIR}/${tumour_id}/${tumour_id}_recalibrated.bam"
        normal_bam="${BAM_DIR}/${normal_id}/${normal_id}_recalibrated.bam"

        if [ ! -f "${tumour_bam}" ]; then
            echo "Error: Tumour BAM file not found: ${tumour_bam}"
            return 1
        fi
        
        if [ ! -f "${normal_bam}" ]; then
            echo "Error: Normal BAM file not found: ${normal_bam}"
            return 1
        fi

        samtools index "${tumour_bam}"
        samtools index "${normal_bam}"

        ## Run TelomereHunter
        ## * Disable the plotting to prevent R errors
        telomerehunter \
            -ibt "${tumour_bam}" \
            -ibc "${normal_bam}" \
            -p "${tumour_id}" \
            -b "${BANDING_FILE}" \
            -o "${TELOMERHUNTER_DIR}" \
            -pl \
            --noFiltering \
            --plotNone
            # >& "${TELOMERHUNTER_DIR}/${tumour_id}/${tumour_id}.telomerehunter.log"
    fi

}

export -f run_telomerhunter

tumour_ids=$(find "${BAM_DIR}" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | grep "T" | sort)

PARALLEL_JOBS=20

echo "${tumour_ids}" | parallel \
    --jobs "$PARALLEL_JOBS" \
    run_telomerhunter {}