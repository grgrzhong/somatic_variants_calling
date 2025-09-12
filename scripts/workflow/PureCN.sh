#!/bin/bash
#SBATCH --job-name=somatic_variants_calling
#SBATCH --partition=amd
#SBATCH --time=72:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.out
#SBATCH --error=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

#############################################################################
# PureCN Workflow Script
# from: https://bioconductor.org/packages/release/bioc/vignettes/PureCN/inst/doc/Quick.html
# Author: Guorui Zhong 
# Date: 2025-09-08
###############################################################################

# conda activate renv
# export PURECN="/home/zhonggr/miniforge3/envs/renv/lib/R/library/PureCN/extdata"
# Rscript $PURECN/PureCN.R --help
# apptainer shell \
#     ./containers/purecn.sif \
#     Rscript "$PURECN/PureCN.R" --help

export REFERENCE_DIR="/lustre1/g/path_my/Reference"
export PROJECT_DIR="/lustre1/g/path_my/pipeline/somatic_variants_calling"
export CONTAINER_DIR="${PROJECT_DIR}/containers"
export FASTA="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"

export PURECN_REF_DIR="${REFERENCE_DIR}/PureCN"
export PURE_PATH="/usr/local/lib/R/site-library/PureCN/extdata"

export PON_DIR="${REFERENCE_DIR}/PON-Mutect"
export PON="${PON_DIR}/pon.vcf.gz"

export BAM_DIR="${PROJECT_DIR}/data/DFSP/BAM"
export MUTECT2_DIR="${PROJECT_DIR}/data/DFSP/Mutect2"
export PURECN_DIR="${PROJECT_DIR}/data/DFSP/PureCN"

## Jobs to run in parallel
export PARALLEL_JOBS=20

## "=========================================================================="
## Step1: Generate reference files for PureCN
## "=========================================================================="
# echo "$(date +"%F") $(date +"%T") - Generating reference files ..."
# singularity exec \
#     --bind ${PROJECT_DIR} \
#     --bind ${REFERENCE_DIR} \
#     ${CONTAINER_DIR}/purecn.sif \
#     Rscript "${PURE_PATH}/IntervalFile.R" \
#         --in-file ${INTERVAL} \
#         --fasta ${FASTA} \
#         --off-target \
#         --genome hg38 \
#         --force \
#         --out-file $PURECN_REF_DIR/baits_hg38_intervals.txt \
#         --export $PURECN_REF_DIR/baits_optimized_hg38.bed \
#         >& ${PURECN_REF_DIR}/IntervalFile.log

## "=========================================================================="
## Step2: Calculate GC-normalized coverages:
## "=========================================================================="
## For each sample, tumor and normal, calculate GC-normalized coverages:
# export sample_id=DFSP-001-T
# mkdir -p ${PURECN_DIR}/${sample_id}
echo "$(date +"%F") $(date +"%T") - Calculating GC-normalized coverages ..."

calculate_coverage() {
    
    local sample_id=$1
    
    ## Create output directory if it doesn't exist
    mkdir -p "${PURECN_DIR}/${sample_id}"

    echo "$(date +"%F") $(date +"%T") - Processing ${sample_id} ..."
    
    ## Calculate GC-normalized coverages
    singularity exec \
        --bind ${PROJECT_DIR} \
        --bind ${REFERENCE_DIR} \
        ${CONTAINER_DIR}/purecn.sif \
        Rscript "${PURE_PATH}/Coverage.R" \
            --bam "${BAM_DIR}/${sample_id}/${sample_id}_recalibrated.bam" \
            --out-dir "${PURECN_DIR}/${sample_id}" \
            --intervals "${PURECN_REF_DIR}/baits_hg38_intervals.txt" \
            --cores 4 \
            >& "${PURECN_DIR}/${sample_id}/coverage.log"
}

export -f calculate_coverage

sample_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)

echo "${sample_ids}" | parallel \
    --jobs "$PARALLEL_JOBS" \
    calculate_coverage {}

## "=========================================================================="
## Step3: Build a normal database for coverage normalization:
## "=========================================================================="
# ls -a $PURECN_DIR/normals/*_coverage.txt.gz | cat > ${PURECN_REF_DIR}/normal_coverages.list
# Create a list of normal coverage files

echo "$(date +"%F") $(date +"%T") - Building normal database ..."

find "${PURECN_DIR}" -name "*-N_coverage.txt.gz" > "${PURECN_REF_DIR}/normal_coverages.list"

singularity exec \
    --bind ${PROJECT_DIR} \
    --bind ${REFERENCE_DIR} \
    ${CONTAINER_DIR}/purecn.sif \
    Rscript "$PURE_PATH/NormalDB.R" \
        --out-dir "${PURECN_REF_DIR}" \
        --coverage-files "${PURECN_REF_DIR}/normal_coverages.list" \
        --normal-panel "${PON}" \
        --genome hg38

## "=========================================================================="
## Step4: Run the main PureCN analysis for each tumor sample
## "=========================================================================="
echo "$(date +"%F") $(date +"%T") - PureCN analysis ..."
purecn_analysis() {

    local tumor_id=$1

    echo "$(date +"%F") $(date +"%T") - Processing ${tumor_id} ..."

    singularity exec \
        --bind ${PROJECT_DIR} \
        --bind ${REFERENCE_DIR} \
        ${CONTAINER_DIR}/purecn.sif \
        Rscript "$PURE_PATH/PureCN.R" \
            --out "${PURECN_DIR}/${tumor_id}" \
            --tumor "${PURECN_DIR}/${tumor_id}/${tumor_id}_coverage_loess.txt.gz" \
            --sampleid "${tumor_id}" \
            --vcf "${MUTECT2_DIR}/${tumor_id}/${tumor_id}.mutect2.vcf.gz" \
            --normaldb ${PURECN_REF_DIR}/normalDB_hg38.rds \
            --intervals ${PURECN_REF_DIR}/baits_hg38_intervals.txt \
            --genome hg38

}

tumour_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep "T" | sort)

echo "${tumour_ids}" | parallel \
    --jobs "$PARALLEL_JOBS" \
    purecn_analysis {}