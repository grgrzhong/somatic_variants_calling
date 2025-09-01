#!/bin/bash

#############################################################################
# This script runs GISTIC2 for copy number analysis on WES data.
# To find out the recurrence of copy number alterations in a set of samples.
###############################################################################

## Setup conda env
# conda create -n gistic2 -y
# conda install hcc::gistic2 
source "$(dirname "${BASH_SOURCE[0]}")/conf/config.sh"

# Activate the conda environment
conda activate gistic2

export GISTIC2_DIR="/mnt/f/projects/250224_sarcoma_multiomics/data/wes/GISTIC2/somatic_matched"
export GISTIC_QVL_THRESHOLD=0.25
export GISTIC_PARALLEL_JOBS=6

## create output directory if it does not exist
mkdir -p "$GISTIC2_DIR"

## Test file
# segment_file="/mnt/f/projects/sarcoma_multiomics/data/wes/GISTIC2/FS-DFSP/FS-DFSP.tsv"
# segment_file="/mnt/f/projects/250224_sarcoma_multiomics/data/temp/GISTIC2/FGFR_mskimpact_segments/FGFR_mskimpact_segments.tsv"
# segment_file="/mnt/f/projects/250224_sarcoma_multiomics/data/wes/GISTIC2/somatic_matched/all_tumors/all_tumors.tsv"
# segment_file="/mnt/f/projects/250224_sarcoma_multiomics/data/wes/test/all_tumors/all_tumors2.tsv"
# segment_file="/mnt/f/projects/250224_sarcoma_multiomics/data/wes/test/somatic_matched/FS-DFSP/FS-DFSP.tsv"
# segment_file="/mnt/f/projects/250224_sarcoma_multiomics/data/wes/GISTIC2/somatic_matched/FST/FST.tsv"

# Function to run GISTIC2 on one file
run_gistic2() {

    local segment_file=$1
    
    # Extract the base name of the segment file
    group=$(basename "${segment_file}" .tsv)

    echo "$(date +"%F") $(date +"%T") - (${group}) Running GISTIC2 with q-value ${GISTIC_QVL_THRESHOLD} ..."
    
    # Create output directory for GISTIC2 results
    output_dir="${GISTIC2_DIR}/${group}/qval_${GISTIC_QVL_THRESHOLD}"
    mkdir -p "$output_dir"

    # Run GISTIC2
    # 0.25 - default threshold for q-value
    # 0.1  - Allows more lenient detection of significant peaks
    # 0.05 - More conservative, fewer false positives

    gistic2 \
        -b "${output_dir}" \
        -seg "${segment_file}" \
        -refgene "${GISTIC2_REFERENCE}" \
        -qvt "${GISTIC_QVL_THRESHOLD}" \
        -ta 0.3 \
        -td 0.3 \
        -conf 0.99 \
        >& "${output_dir}/gistic2.log"
        # -ta 0.1 \
        # -td 0.1 \
        # -armpeel 1 \
        # -js 8 \
        # -genegistic 1 \
        # -conf 0.99 \
        # -brlen 0.8 \
}

# Export the function for GNU Parallel
export -f run_gistic2

# Find all ready-to-use segment files and run in parallel
find "$GISTIC2_DIR" -name "*.tsv" | 
    parallel \
    --jobs ${GISTIC_PARALLEL_JOBS} \
    run_gistic2 {}
