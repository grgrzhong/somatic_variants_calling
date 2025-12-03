#!/bin/bash
#SBATCH --job-name=Strelka2
#SBATCH --partition=amd
#SBATCH --time=6:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.out
#SBATCH --error=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

## Tested only works with this specific build of strelka
# conda create -n strelka strelka=2.9.10=h9ee0642_1

source $(conda info --base)/etc/profile.d/conda.sh
conda activate strelka

## Inputs and outputs
export PROJECT_DIR="/lustre1/g/path_my/pipeline/somatic_variants_calling"
export CONTAINER_DIR="${PROJECT_DIR}/containers"
export REFERENCE_DIR="/lustre1/g/path_my/Reference"
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
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
    ## Extract patient_id from tumour_id (everything before the second hyphen)
    patient_id=$(echo "$tumour_id" | cut -d'-' -f1,2)
    
    ## Use matched normal or defined normal sample
    normal_id="${patient_id}-N"
    
    ## Temporary output directory (WSL can not create makeinfo)

    # Check if presence of paired normal samples
    if [ -d "${BAM_DIR}/${normal_id}" ]; then
        
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

        ## Run Strelka configuration
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Strelka configuration ..."
        
        configureStrelkaSomaticWorkflow.py \
            --normalBam "${normal_bam}" \
            --tumorBam "${tumour_bam}" \
            --callRegions "${INTERVAL}" \
            --exome \
            --indelCandidates "${SV_CALL_DIR}/${tumour_id}/Manta/${tumour_id}.candidateSmallIndels.vcf.gz" \
            --referenceFasta "${REFERENCE}" \
            --runDir "${out_dir}" \
            >& "${out_dir}/strelka_configure.log"

        ## Run the Strelka workflow
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Strelka workflow ..."
        "${out_dir}/runWorkflow.py" \
            -m local \
            -j 16 \
            >& "${out_dir}/strelka_workflow.log"

        ## Combine the SNV and Indel VCFs
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Combining VCFs ..."
        bcftools concat -Oz -a \
            "${out_dir}/results/variants/somatic.snvs.vcf.gz" \
            "${out_dir}/results/variants/somatic.indels.vcf.gz" \
            -o "${out_dir}/results/variants/${tumour_id}_strelka.vcf.gz"
        
        ## Normalize the combined VCF
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Normalizing VCF ..."
        bcftools norm -m-both -f "${REFERENCE}" -Oz \
            -o "${out_dir}/results/variants/${tumour_id}_strelka_normalized.vcf.gz" \
            "${out_dir}/results/variants/${tumour_id}_strelka.vcf.gz"
        tabix "${out_dir}/results/variants/${tumour_id}_strelka_normalized.vcf.gz"

    fi
}

export -f strelka

## "==========================================================================="
## Run Strelka for each tumour sample
## "==========================================================================="
# tumour_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep "T" | sort)
tumour_ids="SARC-035-T"
PARALLEL_JOBS=1

echo "${tumour_ids}" | parallel --jobs "$PARALLEL_JOBS" strelka {}