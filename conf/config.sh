#!/bin/bash

#############################################################################
# Somatic Mutation Workflow - Configuration Script
# This script sets up the environment and variables for the somatic mutation analysis
# workflow.
# It configures directories, container paths, and reference files.
#############################################################################

# Activate conda environment
# source $(conda info --base)/etc/profile.d/conda.sh
# conda activate apptainer

# Set project directory relative to this script
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PROJECT_DIR

# Module directory containing scripts for processing data
export MODULE_DIR="${PROJECT_DIR}/scripts/modules"

## "==========================================================================="
## Input options
## "==========================================================================="
export INPUT_DIR="${1:-${PROJECT_DIR}/data/test_data/RAW}"
export OUTPUT_DIR="${2:-${PROJECT_DIR}/data/test_data}"

## Number of jobs to run in parallel, must be less than the number of samples
export PARALLEL_JOBS="${3:-8}"
export REFERENCE_DIR="${4:-/lustre1/g/path_my/Reference}"
export SOFTWARE_DIR="${5:-/lustre1/g/path_my/Software}"
export CONTAINER_DIR="${6:-${PROJECT_DIR}}/containers"

## "==========================================================================="
## Output directories and Reference files
## "==========================================================================="
mkdir -p "$OUTPUT_DIR"

export FASTQ_TRIM_DIR="${OUTPUT_DIR}/Fastq-trimmed"
export FASTQC_TRIM_DIR="${OUTPUT_DIR}/FastQC-trimmed"

export BAM_DIR="${OUTPUT_DIR}/BAM"

export MUTECT2_DIR="${OUTPUT_DIR}/Mutect2"

export CNV_FACETS_DIR="${OUTPUT_DIR}/CNV_FACETS"

export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export GERMLINE="${REFERENCE_DIR}/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
export PON_DIR="${REFERENCE_DIR}/PON-Mutect"
export PON="${PON_DIR}/pon.vcf.gz"
export PON_PCGR="${PON_DIR}/pon_pcgr.vcf.gz"

export PCGR_DIR="${OUTPUT_DIR}/PCGR"
export PCGR_REFERENCE="${REFERENCE_DIR}/PCGR_reference/20250314"
export VEP_CACHE="${REFERENCE_DIR}/VEP_cache"

export ANNOVAR_DIR="${OUTPUT_DIR}/Annovar"

export FUNOCATOR_ANNOTATION_FILE="${REFERENCE_DIR}/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/"
export ANNOTATION="${REFERENCE_DIR}/Gencode/annotation_protein_coding.bed"

export DBSNP="${REFERENCE_DIR}/Population_database/dbSNP.vcf.gz"
export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
export BAIT_INTERVALS="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-probes-hg38.interval_list"
export TARGET_INTERVALS="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-targets-hg38.interval_list"

export SV_CALL_DIR="${OUTPUT_DIR}/SV-call"

# Print out the environment information
echo "========================================================================"
echo "Clinical RNA Fusion Analysis Workflow - Configuration"
echo "========================================================================"
echo "Project Directory:               $PROJECT_DIR"
echo "Reference Directory:             $REFERENCE_DIR"
echo "Software Directory:              $SOFTWARE_DIR"
echo "Container Directory:             $CONTAINER_DIR"
echo "Input Directory:                 $INPUT_DIR"
echo "Output Directory:                $OUTPUT_DIR"
