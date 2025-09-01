#!/bin/bash

#############################################################################
# Somatic Mutation Workflow - BWA Alignment Script
# This script performs BWA alignment, UMI tagging, sorting, marking duplicates,
# base recalibration, and collects alignment statistics.
#############################################################################

# Create the output directories if they do not exist
mkdir -p "${BAM_DIR}"

# Function to process a single sample through fastp and BWA
bwa_mapping() {

    local sample=$1

    fastq_1="${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_1.fastq.gz"
    fastq_2="${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_2.fastq.gz"

    # Check if input files exist
    if [[ ! -f "${fastq_1}" || ! -f "${fastq_2}" ]]; then
        echo "ERROR: Fastq files for ${sample} not found!"
        return 1
    fi

    # Create the output directory for the sample
    mkdir -p "${BAM_DIR}/${sample}"

    # BWA Alignment
    echo "$(date +"%F") $(date +"%T") - (${sample}) Aligning to reference genome ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bwa.sif" \
        bwa mem \
        -M \
        -t 8 \
        -R "@RG\tID:${sample}\tLB:XGenV2\tPL:ILLUMINA\tPM:NOVASEQ\tSM:${sample}\tPU:NA" \
        "${REFERENCE}" \
        "${fastq_1}" \
        "${fastq_2}" >"${BAM_DIR}/${sample}/${sample}.bwa.sam"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: BWA alignment failed for ${sample}"
        return 1
    fi

    # Convert SAM to BAM
    echo "$(date +"%F") $(date +"%T") - (${sample}) Converting SAM to BAM ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/samtools.sif" \
        samtools \
        view \
        -Sb "${BAM_DIR}/${sample}/${sample}.bwa.sam" >"${BAM_DIR}/${sample}/${sample}.bwa.bam"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: SAM to BAM conversion failed for ${sample}"
        return 1
    fi

    # Extract and tag UMI
    echo "$(date +"%F") $(date +"%T") - (${sample}) Extract and tag UMI ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/pysam.sif" \
        python "${MODULE_DIR}/tag_umi.py" \
        -I "${BAM_DIR}/${sample}/${sample}.bwa.bam" \
        -O "${BAM_DIR}/${sample}/${sample}.umi.bam"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: UMI tagging failed for ${sample}"
        return 1
    fi

    # Sort by coordinate
    echo "$(date +"%F") $(date +"%T") - (${sample}) Sorting BAM file ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/samtools.sif" \
        samtools sort \
            "${BAM_DIR}/${sample}/${sample}.umi.bam" \
            -o "${BAM_DIR}/${sample}/${sample}.sorted.bam"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: BAM sorting failed for ${sample}"
        return 1
    fi

    # Mark duplicates
    echo "$(date +"%F") $(date +"%T") - (${sample}) Marking duplicates ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk.sif" \
        gatk --java-options "-Xmx4g" MarkDuplicates \
        -I "${BAM_DIR}/${sample}/${sample}.sorted.bam" \
        -M "${BAM_DIR}/${sample}/${sample}.metrics.txt" \
        -O "${BAM_DIR}/${sample}/${sample}.marked.bam" \
        --BARCODE_TAG "RX"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Marking duplicates failed for ${sample}"
        return 1
    fi

    # Index BAM
    echo "$(date +"%F") $(date +"%T") - (${sample}) Indexing BAM file ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/samtools.sif" \
        samtools index "${BAM_DIR}/${sample}/${sample}.marked.bam"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: BAM indexing failed for ${sample}"
        return 1
    fi

    # Base recalibration
    echo "$(date +"%F") $(date +"%T") - (${sample}) Running base recalibration ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk.sif" \
        gatk BaseRecalibrator \
        -I "${BAM_DIR}/${sample}/${sample}.marked.bam" \
        -R "${REFERENCE}" \
        -L "${INTERVAL}" \
        -O "${BAM_DIR}/${sample}/${sample}.recal.table" \
        --known-sites "${DBSNP}"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Base recalibration failed for ${sample}"
        return 1
    fi

    # Apply BQSR
    echo "$(date +"%F") $(date +"%T") - (${sample}) Applying BQSR ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk.sif" \
        gatk ApplyBQSR \
            -I "${BAM_DIR}/${sample}/${sample}.marked.bam" \
            -O "${BAM_DIR}/${sample}/${sample}.recal.bam" \
            -L "${INTERVAL}" \
            -bqsr "${BAM_DIR}/${sample}/${sample}.recal.table" \
            --create-output-bam-md5

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Applying BQSR failed for ${sample}"
        return 1
    fi

    # Collect HS metrics
    echo "$(date +"%F") $(date +"%T") - (${sample}) Collecting HsMetrics ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk.sif" \
        gatk CollectHsMetrics \
            -I "${BAM_DIR}/${sample}/${sample}.recal.bam" \
            -O "${BAM_DIR}/${sample}/${sample}.hsmetrics.txt" \
            -R "${REFERENCE}" \
            -BI "${BAIT_INTERVAL}" \
            -TI "${TARGET_INTERVAL}"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Collecting HsMetrics failed for ${sample}"
        return 1
    fi

    # Generate alignment stats
    echo "$(date +"%F") $(date +"%T") - (${sample}) Generating alignment stats ..."
    
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bamtools.sif" \
        bamtools stats \
        -in "${BAM_DIR}/${sample}/${sample}.recal.bam" \
        -out "${BAM_DIR}/${sample}/${sample}.stats.txt"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Generating alignment stats failed for ${sample}"
        return 1
    fi

    # Clean up intermediate files
    echo "$(date +"%F") $(date +"%T") - (${sample}) Cleaning up intermediate files ..."
    rm -f  \
        "${BAM_DIR}/${sample}/${sample}.bwa.sam" \
        "${BAM_DIR}/${sample}/${sample}.bwa.bam" \
        "${BAM_DIR}/${sample}/${sample}.umi.bam" \
        "${BAM_DIR}/${sample}/${sample}.sorted.bam" \
        "${BAM_DIR}/${sample}/${sample}.marked.bam" \
        "${BAM_DIR}/${sample}/${sample}.marked.bam.bai"
}

# Export function to make it available to GNU parallel
export -f bwa_mapping

# Get unique sample names from trimmed fastq files
samples=$(find "${FASTQ_TRIM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n")

# Process samples in parallel
echo "$samples" |
    parallel \
        --jobs "$PARALLEL_JOBS" \
        bwa_mapping {}
