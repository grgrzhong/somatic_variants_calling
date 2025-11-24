#!/bin/bash
#SBATCH --job-name=Somatic_Call
#SBATCH --partition=amd
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.out
#SBATCH --error=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

###############################################################################
## Authors: Zhong Guorui
## Date: 2025-06-10
## Key features:
##      1. Added parallelization to run each steps for multiple samples;
##      2. Use singularity container to keep the environment consistent;
##      3. Only need to install parallel and apptainer
###############################################################################

## Determine the script directory - works both locally and on SLURM
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    ## Running on SLURM - use the submit directory
    PIPELINE_DIR="${SLURM_SUBMIT_DIR}"
else
    ## Running locally - use the directory containing this script
    PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

## Load configuration
source "${PIPELINE_DIR}/conf/config.sh" "$@"

## Load configuration
# source "$(dirname "${BASH_SOURCE[0]}")/conf/config.sh"

## Add better error handling
set -e  # Exit on any error
set -u  # Exit on undefined variables
set -o pipefail  # Exit if any command in pipeline fails

## Log the start time to calculate total time taken
start_time=$(date +"%F %T")

echo "======================================================================="
echo "Somatic variants calling Workflow - Pipeline"
echo "======================================================================="

## Step 1: Preprocessing
echo "$(date +"%F") $(date +"%T") Step 1: Preprocessing and QC ..."
bash "${PROJECT_DIR}/scripts/workflow/step_01_preprocessing.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: Preprocessing steps failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 1: Preprocessing and QC (✓) "

## Step 2: BWA alignment
echo "$(date +"%F") $(date +"%T") Step 2: BWA-mem alignment ..."
bash "${PROJECT_DIR}/scripts/workflow/step_02_bwa_alignment.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: BWA alignment steps failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 2: BWA alignment (✓) "

## Step 3: Mutect2 call
echo "$(date +"%F") $(date +"%T") Step 3: Mutect2 calling ..."
bash "${PROJECT_DIR}/scripts/workflow/step_03_mutect2_call.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: Mutect2 calling steps failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 3: Mutect2 calling (✓) "

## Step4: CNV FACETS
echo "$(date +"%F") $(date +"%T") Step 4: CNV FACETS calling ..."
bash "${PROJECT_DIR}/scripts/workflow/step_04_cnv_facets.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: CNV FACETS calling steps failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 4: CNV FACETS calling (✓) "

## Step5: PCGR annotation
echo "$(date +"%F") $(date +"%T") Step 5: PCGR annotation ..."
bash "${PROJECT_DIR}/scripts/workflow/step_05_pcgr_annotation.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: PCGR annotation steps failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 5: PCGR annotation (✓) "

## Step6: Annovar annotation
echo "$(date +"%F") $(date +"%T") Step 6: Annovar annotation ..."
bash "${PROJECT_DIR}/scripts/workflow/step_06_annovar_annotation.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: Annovar annotation steps failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 6: Annovar annotation (✓) "

## Step7: SV calling
echo "$(date +"%F") $(date +"%T") Step 7: SV calling ..."
bash "${PROJECT_DIR}/scripts/workflow/step_07_sv_call.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: SV calling steps failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 7: SV calling (✓) "

## Log the time of completion
end_time=$(date +"%F %T")
echo "$(date +"%F") $(date +"%T") Pipeline completed at: $end_time"

## Calculate total time taken
start_time_sec=$(date -d "$start_time" +%s)
end_time_sec=$(date -d "$end_time" +%s)
total_time=$((end_time_sec - start_time_sec))
hours=$((total_time / 3600))
minutes=$(((total_time % 3600) / 60))
seconds=$((total_time % 60))

## Print the pipeline summary
echo "======================================================================="
echo "Somatic variants calling Workflow - Running Summary"
echo "======================================================================="
echo "Start Time:               $start_time"
echo "End Time:                 $end_time"
echo "Total Time:               ${hours}h ${minutes}m ${seconds}s"
echo "Output directory:         $OUTPUT_DIR"