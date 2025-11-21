#!/bin/bash


source $(conda info --base)/etc/profile.d/conda.sh
# conda create -n medicc
# conda install -c bioconda -c conda-forge medicc2
conda activate medicc

export MEDICC_DIR="/lustre1/g/path_my/pipeline/somatic_variants_calling/data/DFSP/MEDICC2"

case_ids=$(find "${MEDICC_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n")

## Loop through each case
for case_id in $case_ids; do
    
    input_file="${MEDICC_DIR}/${case_id}/${case_id}.medicc_input.tsv"
    output_dir="${MEDICC_DIR}/${case_id}"
    
    ## Check if input file exists
    if [[ ! -f "$input_file" ]]; then
        echo "ERROR: Input file not found for case ${case_id}: $input_file" >&2
        continue
    fi

    ## Run MEDICC2 for the case
    echo "$(date +"%F") $(date +"%T") - (${case_id}) Running MEDICC2 ..."
    
    medicc2 "${input_file}" "${output_dir}" --prefix "${case_id}" --events \
        >& "${output_dir}/${case_id}.medicc2.log"

done
