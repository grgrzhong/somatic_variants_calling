#!/bin/bash

#############################################################################
## Somatic Mutation Workflow - Mutect2 calling
## Author:  Zhong Guorui
#############################################################################

## Create the output directories if they do not exist
mkdir -p "${MUTECT2_DIR}"

## "==========================================================================="
## Function to run mutect2 and filtering
## "==========================================================================="
mutect2_call() {
    
    local tumour_id=$1
    
    # Extract patient_id from tumour_id (everything before the second hyphen)
    patient_id=$(echo "${tumour_id}" | cut -d'-' -f1,2)
    
    normal_id=${patient_id}-N

    if [[ -f "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.final.vcf.gz" ]]; then
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Mutect2 calling already done. Skipping ..."
        return 0
    fi

    # Create output directory
    mkdir -p "${MUTECT2_DIR}/${tumour_id}"

    ## ========================================================================
    ## Step1. Get Pileup Summaries
    ## ========================================================================
    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Getting Pileup Summaries ..."
    
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk4.sif" \
        gatk GetPileupSummaries \
            -I "${BAM_DIR}/${tumour_id}/${tumour_id}_recalibrated.bam" \
            -V "${REFERENCE_DIR}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz" \
            -L "${REFERENCE_DIR}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz" \
            -O "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.getpileupsummaries.table" \
            >& "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.getpileupsummaries.log"

    # Get Pileup Summaries if presence of paired normal samples
    if [ -d "${BAM_DIR}/${normal_id}" ]; then
        echo "$(date +"%F") $(date +"%T") - (${normal_id}) Getting Pileup Summaries ..."
    
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk4.sif" \
        gatk GetPileupSummaries \
            -I "${BAM_DIR}/${normal_id}/${normal_id}_recalibrated.bam" \
            -V "${REFERENCE_DIR}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz" \
            -L "${REFERENCE_DIR}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz" \
            -O "${MUTECT2_DIR}/${tumour_id}/${normal_id}.getpileupsummaries.table" \
            >& "${MUTECT2_DIR}/${tumour_id}/${normal_id}.getpileupsummaries.log"

    fi

    ## ====================================================================
    ## Step2. Calculate contamination
    ## ====================================================================
    ## Calculate contamination based on tumour samples only
    if [ -d "${BAM_DIR}/${normal_id}" ]; then

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Calculating contamination of matched samples ..."

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/gatk4.sif" \
            gatk CalculateContamination \
                -I "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.getpileupsummaries.table" \
                -matched "${MUTECT2_DIR}/${tumour_id}/${normal_id}.getpileupsummaries.table" \
                -O "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.contamination.table" \
                -segments "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.segments.table" \
                >& "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.contamination.log"

    else
        # Calculate contamination based on tumour samples only
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Calculating contamination of unmatched samples ..."
        
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/gatk4.sif" \
            gatk CalculateContamination \
                -I "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.getpileupsummaries.table" \
                -O "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.contamination.table" \
                -segments "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.segments.table" \
                >& "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.contamination.log"

    fi

    ## ====================================================================
    ## Step3. Mutect2 call variants
    ## ====================================================================
    # Call variant on paired tumour samples
    if [ -d "${BAM_DIR}/${normal_id}" ]; then

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Mutect2 call on matched samples ..."

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/gatk4.sif" \
            gatk --java-options -Xmx4g Mutect2 \
                -I "$BAM_DIR/${tumour_id}/${tumour_id}_recalibrated.bam" \
                -I "$BAM_DIR/${normal_id}/${normal_id}_recalibrated.bam" \
                -normal "$normal_id" \
                -R "${REFERENCE}" \
                -L "${INTERVAL}" \
                --germline-resource "${GERMLINE}" \
                --panel-of-normals "${PON}" \
                --f1r2-tar-gz "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.f1r2.tar.gz" \
                --native-pair-hmm-threads 8 \
                --callable-depth 20 \
                -O "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.mutect2.vcf.gz" \
                -bamout "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.realigned.bam" \
                >& "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.mutect2.log"

    else
        ## Mutect2 call variant on unpaired tumour samples
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Mutect2 call on unmatched samples ..."

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/gatk4.sif" \
            gatk --java-options -Xmx8g Mutect2 \
                -I "${BAM_DIR}/${tumour_id}/${tumour_id}_recalibrated.bam" \
                -R "${REFERENCE}" \
                -L "${INTERVAL}" \
                --germline-resource "${GERMLINE}" \
                --panel-of-normals "${PON}" \
                --f1r2-tar-gz "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.f1r2.tar.gz" \
                --callable-depth 20 \
                -O "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.mutect2.vcf.gz" \
                -bamout "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.realigned.bam" \
                >& "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.mutect2.log"

    fi

    ## ====================================================================
    ## Step4. Learn Read Orientation Model bias for artifacts
    ## ====================================================================
    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Learning Read Orientation Model ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk4.sif" \
        gatk LearnReadOrientationModel \
            -I "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.f1r2.tar.gz" \
            -O "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.readorientationmodel.tar.gz" \
            >& "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.learnreadorientationmodel.log"

    ## ====================================================================
    ## Step5. Learn Read Orientation Model bias for artifacts
    ## ====================================================================
    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Filtering mutect2 calls ..."
    
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk4.sif" \
        gatk --java-options -Xmx4g FilterMutectCalls \
            --variant "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.mutect2.vcf.gz" \
            --stats "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.mutect2.vcf.gz.stats" \
            --reference "${REFERENCE}" \
            --ob-priors "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.readorientationmodel.tar.gz" \
            --contamination-table "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.contamination.table" \
            --tumor-segmentation "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.segments.table" \
            --min-allele-fraction 0.01 \
            --unique-alt-read-count 1 \
            --output "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.filtermutectcalls.vcf.gz" \
            >& "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.filtermutectcalls.log"

    ## ====================================================================
    ## Step6. Normalize the variants
    ## ====================================================================
    ## Standardize the vcf, split multiallelic variants, and left align indels
    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Normalizing mutect2 calls ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bcftools.sif" \
        bcftools norm \
            "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.filtermutectcalls.vcf.gz" \
            -m-both -f "${REFERENCE}" \
            -Oz \
            -o "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.normalized.vcf.gz"

    ## ====================================================================
    ## Step7. Filter out low-quality or failed variant calls
    ## ====================================================================
    ## Keep variants that have passed all filters (low-quality or failed variant calls)
    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Filtering PASS variants ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bcftools.sif" \
        bcftools view \
            -f PASS "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.normalized.vcf.gz" \
            -o "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.passed.vcf.gz"

    ## ====================================================================
    ## Step8. Filter out blacklist and repeatmasker regions
    ## ====================================================================
    ## Annotate repeatmasker and blacklist regions
    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Annotating repeatmasker regions ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bcftools.sif" \
        bcftools annotate \
            "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.passed.vcf.gz" \
            --header-lines "${MUTECT2_DIR}/vcf.rm.header" \
            --annotations "${REFERENCE_DIR}/RepeatMasker.bed.gz" \
            --columns CHROM,FROM,TO,RepeatMasker \
            --output "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.repeatmasker.vcf.gz"

    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Annotating blacklist regions ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bcftools.sif" \
        bcftools annotate \
            "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.repeatmasker.vcf.gz" \
            --header-lines "${MUTECT2_DIR}/vcf.map.header" \
            --annotations "${REFERENCE_DIR}/blacklist.bed.gz" \
            --columns CHROM,FROM,TO,EncodeDacMapability \
        --output-type z \
        --output "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.blacklist.vcf.gz"

    ## Filter out variants in RepeatMasker or Mapability
    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Filtering RepeatMasker and blacklist regions ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bcftools.sif" \
        bcftools filter \
            "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.blacklist.vcf.gz" \
            -e 'INFO/RepeatMasker != "." || INFO/EncodeDacMapability != "."' \
            -Oz \
            -o "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.filtered.vcf.gz"

    ## Remove the header lines from the final VCF file
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bcftools.sif" \
        bcftools annotate \
            --remove INFO/RepeatMasker,INFO/EncodeDacMapability \
            --output-type z \
            --output "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.final.vcf.gz" \
            "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.filtered.vcf.gz"

    # Index the final VCF file
    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Indexing final VCF file ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/tabix.sif" \
        tabix -p vcf "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.final.vcf.gz"

    # Clean up intermediate files
    echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Cleaning up intermediate files ..."
    rm -f "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.normalized.vcf.gz"
    rm -f "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.passed.vcf.gz"
    rm -f "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.repeatmasker.vcf.gz"
    rm -f "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.blacklist.vcf.gz"
    rm -f "${MUTECT2_DIR}/${tumour_id}/${tumour_id}.filtered.vcf.gz"

}

export -f mutect2_call

## "==========================================================================="
## Create header files needed for repeatmasker and blacklist annotation
## "==========================================================================="
echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > "${MUTECT2_DIR}/vcf.rm.header"
echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" > "${MUTECT2_DIR}/vcf.map.header"

## "==========================================================================="
## Run mutect2 calling for each tumour sample
## "==========================================================================="
tumour_samples=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -name "*T*" -print0 | xargs -0 -n 1 basename)

echo "${tumour_samples}" | parallel --jobs "${PARALLEL_JOBS}" mutect2_call {}

# Clean up temporary header files
rm "${MUTECT2_DIR}/vcf.rm.header"
rm "${MUTECT2_DIR}/vcf.map.header"

