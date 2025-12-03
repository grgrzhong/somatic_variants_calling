#!/bin/bash

## Create the output directories
export PROJECT_DIR="/lustre1/g/path_my/pipeline/somatic_variants_calling"
export CONTAINER_DIR="${PROJECT_DIR}/containers"
export REFERENCE_DIR="/lustre1/g/path_my/Reference"
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export SOFTWARE_DIR="/lustre1/g/path_my/Software"
export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed.gz"
export BAM_DIR="${PROJECT_DIR}/data/SARC/BAM"
export SV_CALL_DIR="${PROJECT_DIR}/data/SARC/SV-call"
export STRELKA_DIR="${PROJECT_DIR}/data/SARC/Strelka2"
mkdir -p "${STRELKA_DIR}"

## "==========================================================================="
## Function to run strelka for a sample
## "==========================================================================="
strelka() {
    
    local tumour_id="$1"
    tumour_id="SARC-035-T"
    ## Extract patient_id from tumour_id (everything before the second hyphen)
    patient_id=$(echo "$tumour_id" | cut -d'-' -f1,2)
    
    ## Use matched normal or defined normal sample
    normal_id="${patient_id}-N"
    
    ## Temporary output directory (WSL can not create makeinfo)

    # Check if presence of paired normal samples
    if [ -d "${BAM_DIR}/${normal_id}" ]; then
        
        if [ -f "${STRELKA_DIR}/${tumour_id}/${tumour_id}.vcf.gz" ]; then
            echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Strelka already done. Skipping ..."
            return 0
        fi
        
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Running Strelka ..."
        
        ## "-----------------------------------------------------------------"
        ## Run Strelka
        ## "-----------------------------------------------------------------"
        out_dir="${STRELKA_DIR}/${tumour_id}"
        rm -rf "${out_dir}"
        mkdir -p "${out_dir}"
        
        ## Check the bam files for tumour and normal samples
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

        ## Run Strelka
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${STRELKA_DIR}:${STRELKA_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            ${CONTAINER_DIR}/strelka.sif \
            configureStrelkaSomaticWorkflow.py \
                --normalBam "${normal_bam}" \
                --tumorBam "${tumour_bam}" \
                --callRegions "${INTERVAL}" \
                --exome \
                --referenceFasta "${REFERENCE}" \
                --runDir "${out_dir}" \
                >& "${out_dir}/strelka_configure.log"
                # --indelCandidates "${SV_CALL_DIR}/${tumour_id}/Manta/${tumour_id}.candidateSmallIndels.vcf.gz" \

        ## Run the workflow
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${STRELKA_DIR}:${STRELKA_DIR}" \
            ${CONTAINER_DIR}/strelka.sif \
            "${out_dir}/runWorkflow.py" -m sge -j 1 -g 12 \
                >& "${out_dir}/strelka_workflow.log"

    fi
}

export -f strelka

## "==========================================================================="
## Run Strelka for each tumour sample
## "==========================================================================="
tumour_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep "T" | sort)

echo "${tumour_ids}" | parallel --jobs "$PARALLEL_JOBS" strelka {}