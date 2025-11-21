#!/bin/bash
#SBATCH --job-name=FACETS2N
#SBATCH --partition=amd
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.out
#SBATCH --error=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Determine the script directory - works both locally and on SLURM
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # Running on SLURM - use the submit directory
    PIPELINE_DIR="${SLURM_SUBMIT_DIR}"
else
    # Running locally - use the directory containing this script
    PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Load configuration
source "${PIPELINE_DIR}/conf/config.sh" "$@"

## Create the environment
# conda create -n facets bioconda::r-facets bioconda::snp-pileup r-devtools htslib=1.3.1

# export DBSNP="${REFERENCE_DIR}/Population_database/dbSNP.vcf.gz"
# export DBSNP="${REFERENCE_DIR}/FACETS/facets_reference_snps_hg38_uniform_1kb.chr_style.vcf.gz"
export DBSNP="${REFERENCE_DIR}/Population_database/00-common_all.vcf.gz"
export BAM_DIR="${PROJECT_DIR}/data/DFSP/BAM"
export FACETS2N_DIR="${PROJECT_DIR}/data/DFSP/FACETS2N"
mkdir -p "${FACETS2N_DIR}"

## "==========================================================================="
## Generate the reference snp-pileup
## "==========================================================================="
reference_normals_dir="${FACETS2N_DIR}/reference_normals"
mkdir -p "${reference_normals_dir}"
export SNP_PILEUP="${HOME}/miniforge3/envs/facets/bin/snp-pileup"

# export snp_pileup_path="${HOME}/miniforge3/envs/snp-pileup/bin/snp-pileup"
normal_bams=$(find "${BAM_DIR}" -type f -name "*-N*_recalibrated.bam" | tr '\n' ' ')

# Check if snp-pileup exists
if [ ! -f "${snp_pileup_path}" ]; then
    echo "Error: snp-pileup not found at ${snp_pileup_path}"
    exit 1
fi

# Check if any normal BAMs were found
if [ -z "$normal_bams" ]; then
    echo "Error: No normal BAM files found"
    exit 1
fi

export PATH="${HOME}/htslib/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/htslib-1.22.1/lib:${LD_LIBRARY_PATH}"
export C_INCLUDE_PATH="${HOME}/htslib-1.22.1/include:${C_INCLUDE_PATH}"
export CPLUS_INCLUDE_PATH="${HOME}/htslib-1.22.1/include:${CPLUS_INCLUDE_PATH}"

    # --snp-pileup-path "${snp_pileup_path}" \
Rscript ./scripts/facets2n/inst/extcode/snp-pileup-wrapper.R \
    --output-prefix "${reference_normals_dir}/reference_normals" \
    --vcf-file "${DBSNP}" \
    --unmatched-normal-BAMS "${normal_bams}"


## "==========================================================================="
## Generate the reference snp-pileup
## "==========================================================================="
run_facets() {

    local tumour_id="$1"

    ## Extract patient_id from tumour_id (everything before the second hyphen)
    patient_id=$(echo "$tumour_id" | cut -d'-' -f1,2)
    
    ## Use matched normal or defined normal sample
    normal_id=${patient_id}-N

    if [ -d "${BAM_DIR}/${normal_id}" ]; then

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

        ## Create output directory for the sample
        sample_facets_dir="${FACETS_DIR}/${tumour_id}"
        rm -rf "${sample_facets_dir}"
        mkdir -p "${sample_facets_dir}"

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) - snp-pileup ..."
        ## Run FACETS snp-pileup
        singularity run \
                --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
                --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
                --bind "${FACETS_DIR}:${FACETS_DIR}" \
                --bind "${BAM_DIR}:${BAM_DIR}" \
                --bind /tmp:/tmp \
                "${CONTAINER_DIR}/facets-suite-dev.img" \
                snp-pileup-wrapper.R \
                    --vcf-file "${DBSNP}" \
                    --normal-bam "${normal_bam}" \
                    --tumor-bam "${tumour_bam}" \
                    --verbose \
                    --output-prefix "${sample_facets_dir}/${tumour_id}" \
                    >& "${sample_facets_dir}/${tumour_id}.snp_pileup.log"

        ## FACETS wrapper
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) - run-facets ..."
        # rm -rf "${sample_facets_dir}/${tumour_id}"
        cd "${sample_facets_dir}" || exit 1

        singularity run \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${FACETS_DIR}:${FACETS_DIR}" \
            --bind "${sample_facets_dir}:${sample_facets_dir}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/facets-suite-dev.img" \
            run-facets-wrapper.R \
                --counts-file "${sample_facets_dir}/${tumour_id}.snp_pileup.gz" \
                --directory . \
                --verbose \
                --genome hg38 \
                --sample-id "${tumour_id}" \
                --cval 150 \
                --purity-cval 25 \
                --everything \
                --normal-depth 15 \
                >& "${sample_facets_dir}/${tumour_id}.run_facets.log"
        
        ## move and merge the all output files to upper directory
        mv "${sample_facets_dir}/${tumour_id}"/* "${sample_facets_dir}/" 2>/dev/null || true
        rmdir "${sample_facets_dir}/${tumour_id}" 2>/dev/null || true

        cd "${PROJECT_DIR}" || exit 1

    fi

}

export -f run_facets

## Run FACETS for each tumour sample
tumour_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep "T" | sort)
PARALLEL_JOBS=30

echo "${tumour_ids}" | parallel --jobs "$PARALLEL_JOBS" run_facets {}