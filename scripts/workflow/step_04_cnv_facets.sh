#!/bin/bash

## "==========================================================================="
## Adapted from cnv_facet analysis (https://github.com/dariober/cnv_facets)
## Authors: Zhong Guorui
## Date: 2025-06-10
## Description: This script runs CNV FACETS analysis for paired samples if available.
## and export files for downstream analysis
## Features:
##      1. Uses GNU Parallel for efficient processing of multiple samples
##      2. Reformat the output segment data
##      3. Modified cnv_facets.R to include gene annotation for neutral segments
## "==========================================================================="

## Create the output directories
mkdir -p "${CNV_FACETS_DIR}"

export LOCAL_CNV_FACETS="${MODULE_DIR}/cnv_facets_modified.R"
# export CNV_FACETS_DIR="${PROJECT_DIR}/data/DFSP/CNV_FACETS_test"

## "==========================================================================="
## Function to run CNV FACETS for a sample
## "==========================================================================="
cnv_facets() {
    
    local tumour_id="$1"

    ## Extract patient_id from tumour_id (everything before the second hyphen)
    patient_id=$(echo "$tumour_id" | cut -d'-' -f1,2)
    
    ## Use matched normal or defined normal sample
    normal_id=${patient_id}-N
    
    ## Temporary output directory (WSL can not create makeinfo)

    # Check if presence of paired normal samples
    if [ -d "${BAM_DIR}/${normal_id}" ]; then
        
        if [ -f "${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.vcf.gz" ]; then
            echo "$(date +"%F") $(date +"%T") - (${tumour_id}) CNV FACETS already done. Skipping ..."
            return 0
        fi
        
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Running CNV FACETS ..."
        
        ## "-----------------------------------------------------------------"
        ## Run CNV FACETS
        ## "-----------------------------------------------------------------"
        out_dir="${CNV_FACETS_DIR}/${tumour_id}"
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

        prefix="${out_dir}/${tumour_id}"

        ## Run FACETS
        ## Note: when runing the singularity, should not use the mounted path 
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind "${LOCAL_CNV_FACETS}:/opt/conda/bin/cnv_facets.R" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/cnv_facets.sif" \
            cnv_facets.R \
                --snp-tumour "${tumour_bam}" \
                --snp-normal "${normal_bam}" \
                --snp-vcf "${DBSNP}" \
                --snp-nprocs 16 \
                --cval 25 150 \
                --targets "${INTERVAL}" \
                --annotation "${ANNOTATION}" \
                --out "${prefix}" \
                >& "${prefix}.facets.log"
        
        ## if it returns error, then skip following steps
        if [ $? -ne 0 ]; then
            echo "$(date +"%F") $(date +"%T") - (${tumour_id}) CNV FACETS failed. Skipping ..."
            echo "${tumour_id}" >> "${CNV_FACETS_DIR}/cnv_facets_failed_samples.txt"
            return 1
        fi

        ## "-----------------------------------------------------------------"
        ## Convert the VCF output to TSV format, Without any filtering
        ## "-----------------------------------------------------------------"
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Converting VCF to TSV ..."
        facets_segment="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.facets_segment.tsv"
        facets_gene="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.facets_gene.tsv"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${CNV_FACETS_DIR}:${CNV_FACETS_DIR}" \
            --bind "${MODULE_DIR}:${MODULE_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/r.sif" \
            Rscript "${MODULE_DIR}/cnv_facets_convert_vcf_to_tsv.R" \
                --input "${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.vcf.gz" \
                --output_segment "${facets_segment}" \
                --output_gene "${facets_gene}"

        ## "-----------------------------------------------------------------"
        ## Prepare the Copy number segments for pgcr input
        ## Note: coordinates must be one-based; 
        ## https://sigven.github.io/pcgr/articles/input.html
        ## "-----------------------------------------------------------------"
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Exporting segments for PCGR input ..."
        
        pcgr_tsv="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.pcgr.tsv"
        
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${CNV_FACETS_DIR}:${CNV_FACETS_DIR}" \
            --bind "${MODULE_DIR}:${MODULE_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/r.sif" \
            Rscript "${MODULE_DIR}/cnv_facets_export_segments_to_pcgr.R" \
                --input "${facets_segment}" \
                --output "${pcgr_tsv}"

        ## "-----------------------------------------------------------------"
        ## Prepare the segments for GISTIC2 input
        ## "-----------------------------------------------------------------"
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Exporting segments for GISTIC2 input ..."
        
        gistic2_tsv="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.gistic2.tsv"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${CNV_FACETS_DIR}:${CNV_FACETS_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/r.sif" \
            Rscript "${MODULE_DIR}/cnv_facets_export_segments_to_gistic2.R" \
                --input "${facets_segment}" \
                --output "${gistic2_tsv}"

    fi
}

export -f cnv_facets

## "==========================================================================="
## Run FACETS for each tumour sample
## "==========================================================================="
## Clear previous failed samples log
rm -f "${CNV_FACETS_DIR}/cnv_facets_failed_samples.txt"

tumour_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep "T" | sort)

echo "${tumour_ids}" | parallel --jobs "$PARALLEL_JOBS" cnv_facets {}