#!/bin/bash

###############################################################################
## Adapted from PCGR annotation reference (https://sigven.github.io/pcgr/index.html)
## Authors: Zhong Guorui
## Date: 2025-06-10
## Key features: 
##      1. Added prallelization to run PCGR annotation for multiple samples; 
##      2. Use singularity container to keep the environment consistent;
###############################################################################

## Create the output directory
mkdir -p "${PCGR_DIR}"
export COPY_GAIN_THRESHOLD=3

## ===========================================================================
## General configurations
## ===========================================================================

## Create PCGR-compatible PON if it doesn't exist
if [ ! -f "${PON_PCGR}" ]; then
    echo " - Creating PCGR-compatible PoN VCF..."
    { 
        bcftools view -h "${PON}" | head -n -1
        echo '##INFO=<ID=PANEL_OF_NORMALS,Number=0,Type=Flag,Description="Overlap with germline call among panel of normals">'
        bcftools view -h "${PON}" | tail -n 1
        bcftools view -H "${PON}" | awk 'BEGIN{OFS="\t"} {
            if ($8 == ".") {
                $8 = "PANEL_OF_NORMALS"
            } else {
                $8 = $8 ";PANEL_OF_NORMALS"
            }
            print $0
        }'
    } | bgzip > "${PON_PCGR}" && tabix -p vcf "${PON_PCGR}"
    echo " - PCGR-compatible PoN created: ${PON_PCGR}"
fi

# bcftools view -h "${PON_PCGR}" | grep -i "^##"
# bcftools view -h "${PON_PCGR}" | head -20
# bcftools view -h "${PON_PCGR}" | grep -i "^#CHROM"
# bcftools view -h "${PON}" | grep -i "^##"

## Function to run PCGR annotation
pcgr_annotation() {
    
    local tumour_id=$1

    ## Find the matched normal sample
    patient_id=$(echo "${tumour_id}" | cut -d'-' -f1,2)

    normal_id="${patient_id}-N"

    ## Check if input VCF exists
    input_vcf="${MUTECT2_DIR}/${tumour_id}/${tumour_id}.final.vcf.gz"

    if [ ! -f "${input_vcf}" ]; then
        echo "Input VCF not found: ${input_vcf}"
        return 1
    fi
    
    input_cna="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.pcgr.tsv"

    if [ ! -f "${input_cna}" ]; then
        echo "Input CNA file not found: ${input_cna}"
        return 1
    fi

    ## Create output directory
    output_dir="${PCGR_DIR}/${tumour_id}"
    mkdir -p "${output_dir}"

    ## Define reformatted VCF path
    reformatted_vcf="${output_dir}/${tumour_id}.reformatted.vcf.gz"

    ## ========================================================================
    ## Tumour-Normal PCGR Annotation
    ## ========================================================================
    if [ -d "${BAM_DIR}/${normal_id}" ]; then

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) PCGR tumour-normal mode ..."

        ## Reformat VCF 
        singularity exec \
            --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${MODULE_DIR}:${MODULE_DIR}" \
            --bind "/tmp:/tmp" \
            "${CONTAINER_DIR}/pysam.sif" \
            python "${MODULE_DIR}/pcgr_reformat_vcf_tumour_normal.py" \
                -I "${input_vcf}" \
                -O "${reformatted_vcf}"

        ## Index it
        rm -f "${reformatted_vcf}.tbi"
        singularity exec \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            --bind "${output_dir}:${output_dir}" \
            "${CONTAINER_DIR}/tabix.sif" \
            tabix -p vcf "${reformatted_vcf}"
        
        ## Run PCGR annotation
        singularity exec \
            --bind "${PCGR_REFERENCE}:${PCGR_REFERENCE}" \
            --bind "${VEP_CACHE}:${VEP_CACHE}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            "${CONTAINER_DIR}/pcgr_2.2.3.singularity.sif" \
            pcgr \
            --input_vcf "${reformatted_vcf}" \
            --input_cna "${input_cna}" \
            --vep_dir "${VEP_CACHE}" \
            --refdata_dir "${PCGR_REFERENCE}" \
            --output_dir "${output_dir}" \
            --genome_assembly grch38 \
            --sample_id "${tumour_id}" \
            --assay WES \
            --n_copy_gain "${COPY_GAIN_THRESHOLD}" \
            --effective_target_size_mb 34 \
            --tumor_dp_tag TDP \
            --tumor_af_tag TAF \
            --control_dp_tag NDP \
            --control_af_tag NAF \
            --tumor_dp_min 20 \
            --tumor_af_min 0.05 \
            --control_dp_min 10 \
            --control_af_max 0.01 \
            --estimate_tmb \
            --tmb_dp_min 20 \
            --tmb_af_min 0.05 \
            --estimate_msi \
            --estimate_signatures \
            --vcf2maf \
            --ignore_noncoding \
            --force_overwrite \
            >& "${output_dir}/pcgr.log"

    ## ========================================================================
    ## Tumour-Only PCGR Annotation
    ## ========================================================================
    else
        
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) PCGR tumour-only mode ..."
        
        ## Reformat VCF 
        singularity exec \
            --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${MODULE_DIR}:${MODULE_DIR}" \
            --bind "/tmp:/tmp" \
            "${CONTAINER_DIR}/pysam.sif" \
            python "${MODULE_DIR}/pcgr_reformat_vcf_tumour_only.py" \
                -I "${input_vcf}" \
                -O "${reformatted_vcf}"

        ## Index it
        rm -f "${reformatted_vcf}.tbi"
        singularity exec \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            --bind "${output_dir}:${output_dir}" \
            "${CONTAINER_DIR}/tabix.sif" \
            tabix -p vcf "${reformatted_vcf}"
        
        ## Run PCGR annotation
        singularity exec \
            --bind "${PCGR_REFERENCE}:${PCGR_REFERENCE}" \
            --bind "${VEP_CACHE}:${VEP_CACHE}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            --bind "${PON_DIR}:${PON_DIR}" \
            --bind "/tmp:/tmp" \
            "${CONTAINER_DIR}/pcgr_2.2.3.singularity.sif" \
            pcgr \
                --input_vcf "${reformatted_vcf}" \
                --input_cna "${input_cna}" \
                --pon_vcf "${PON_PCGR}" \
                --vep_dir "${VEP_CACHE}" \
                --refdata_dir "${PCGR_REFERENCE}" \
                --output_dir "$output_dir" \
                --genome_assembly grch38 \
                --sample_id "${tumour_id}" \
                --assay WES \
                --n_copy_gain "${COPY_GAIN_THRESHOLD}" \
                --effective_target_size_mb 34 \
                --tumor_only \
                --tumor_dp_tag TDP \
                --tumor_af_tag TAF \
                --tumor_dp_min 20 \
                --tumor_af_min 0.05 \
                --estimate_tmb \
                --tmb_dp_min 20 \
                --tmb_af_min 0.05 \
                --estimate_msi \
                --estimate_signatures \
                --vcf2maf \
                --ignore_noncoding \
                --force_overwrite \
                >& "${output_dir}/pcgr.log"

    fi
}

export -f pcgr_annotation

## Run PCGR annotation for each tumour samples
tumour_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep "T" | sort)

export PARALLEL_JOBS=20

echo "${tumour_ids}" | parallel \
    --jobs "$PARALLEL_JOBS" \
    pcgr_annotation {}
    