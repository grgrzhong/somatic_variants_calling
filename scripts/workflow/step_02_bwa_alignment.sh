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
    
    log_dir="${BAM_DIR}/${sample}/logs"
    
    mkdir -p "${log_dir}"
    
    # BWA Alignment
    echo "$(date +"%F") $(date +"%T") - (${sample}) BWA alignment ..."
    
    bwa_log="${log_dir}/${sample}.bwa.log"
    
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bwa.sif" \
        bwa mem \
        -M \
        -t 16 \
        -R "@RG\tID:${sample}\tLB:XGenV2\tPL:ILLUMINA\tPM:NOVASEQ\tSM:${sample}\tPU:NA" \
        "${REFERENCE}" \
        "${fastq_1}" \
        "${fastq_2}" > "${BAM_DIR}/${sample}/${sample}.bwa.sam" \
        2>> "${bwa_log}"

    # Convert SAM to BAM
    echo "$(date +"%F") $(date +"%T") - (${sample}) Converting SAM to BAM ..."

    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/samtools.sif" \
        samtools \
        view \
        -@ 8 \
        -Sb "${BAM_DIR}/${sample}/${sample}.bwa.sam" >"${BAM_DIR}/${sample}/${sample}.bwa.bam"

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

    # Sort by coordinate
    echo "$(date +"%F") $(date +"%T") - (${sample}) Sorting UMI BAM ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/samtools.sif" \
        samtools sort \
            "${BAM_DIR}/${sample}/${sample}.umi.bam" \
            -o "${BAM_DIR}/${sample}/${sample}.sorted.bam" \
            >& "${log_dir}/${sample}.SortByCoordinate.log"

    # Mark duplicates
    echo "$(date +"%F") $(date +"%T") - (${sample}) Marking duplicates ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk4.sif" \
        gatk --java-options "-Xmx4g" MarkDuplicates \
        -I "${BAM_DIR}/${sample}/${sample}.sorted.bam" \
        -M "${BAM_DIR}/${sample}/${sample}_metrics.txt" \
        -O "${BAM_DIR}/${sample}/${sample}.marked.bam" \
        --BARCODE_TAG "RX" \
        >& "${log_dir}/${sample}.MarkDuplicates.log"

    # Index BAM
    echo "$(date +"%F") $(date +"%T") - (${sample}) Indexing BAM file ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/samtools.sif" \
        samtools index "${BAM_DIR}/${sample}/${sample}.marked.bam"

    # Base recalibration
    echo "$(date +"%F") $(date +"%T") - (${sample}) Base recalibration ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk4.sif" \
        gatk BaseRecalibrator \
        -I "${BAM_DIR}/${sample}/${sample}.marked.bam" \
        -R "${REFERENCE}" \
        -L "${INTERVAL}" \
        -O "${BAM_DIR}/${sample}/${sample}_recal_table.table" \
        --known-sites "${DBSNP}" \
        >& "${log_dir}/${sample}.BaseRecalibrator.log"

    # Apply BQSR
    echo "$(date +"%F") $(date +"%T") - (${sample}) Applying BQSR ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk4.sif" \
        gatk ApplyBQSR \
            -I "${BAM_DIR}/${sample}/${sample}.marked.bam" \
            -O "${BAM_DIR}/${sample}/${sample}_recalibrated.bam" \
            -L "${INTERVAL}" \
            -bqsr "${BAM_DIR}/${sample}/${sample}_recal_table.table" \
            --create-output-bam-md5 \
            >& "${log_dir}/${sample}.ApplyBQSR.log"

    # Collect HS metrics
    echo "$(date +"%F") $(date +"%T") - (${sample}) Collecting HsMetrics ..."
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/gatk4.sif" \
        gatk CollectHsMetrics \
            -I "${BAM_DIR}/${sample}/${sample}_recalibrated.bam" \
            -O "${BAM_DIR}/${sample}/${sample}_hs_metrics.txt" \
            -R "${REFERENCE}" \
            -BI "${BAIT_INTERVALS}" \
            -TI "${TARGET_INTERVALS}" \
            >& "${log_dir}/${sample}.CollectHsMetrics.log"

    # Generate alignment stats
    echo "$(date +"%F") $(date +"%T") - (${sample}) Generating alignment stats ..."
    
    singularity exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${BAM_DIR}:${BAM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/bamtools.sif" \
        bamtools stats \
        -in "${BAM_DIR}/${sample}/${sample}_recalibrated.bam" > "${BAM_DIR}/${sample}/${sample}_aln_stat.txt"

    # Clean up intermediate files to save disk space
    echo "$(date +"%F") $(date +"%T") - (${sample}) Cleaning up intermediate files ..."

    rm -rf "${BAM_DIR}/${sample}/${sample}.bwa.sam" 
    rm -rf "${BAM_DIR}/${sample}/${sample}.bwa.bam"
    rm -rf "${BAM_DIR}/${sample}/${sample}.umi.bam"
    rm -rf "${BAM_DIR}/${sample}/${sample}.sorted.bam"
    rm -rf "${BAM_DIR}/${sample}/${sample}.marked.bam"
    rm -rf "${BAM_DIR}/${sample}/${sample}.marked.bam.bai"
    rm -rf "${BAM_DIR}/${sample}/${sample}_recalibrated.bam.md5"
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
