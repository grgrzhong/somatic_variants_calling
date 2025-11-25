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

## Create the environment
# conda create -n facets2n snp-pileup r-devtools r-optparse r-parallel

export PROJECT_DIR="/lustre1/g/path_my/pipeline/somatic_variants_calling"
export MODULE_DIR="${PROJECT_DIR}/scripts/modules"
export REFERENCE_DIR="/lustre1/g/path_my/Reference"
export DBSNP="${REFERENCE_DIR}/Population_database/dbSNP.vcf.gz"

# export DBSNP="${REFERENCE_DIR}/FACETS/facets_reference_snps_hg38_uniform_1kb.chr_style.vcf.gz"
# export DBSNP="${REFERENCE_DIR}/Population_database/00-common_all.vcf.gz"
export BAM_DIR="${PROJECT_DIR}/data/DFSP/BAM"
export FACETS2N_DIR="${PROJECT_DIR}/data/DFSP/FACETS2N"
mkdir -p "${FACETS2N_DIR}"

## "==========================================================================="
## Reference snp-pileup
## "==========================================================================="
export FACETS2N_REFERENCE_DIR="${FACETS2N_DIR}/reference_normals"
mkdir -p "${FACETS2N_REFERENCE_DIR}"

echo "$(date +"%F") $(date +"%T") - Reference normals snp-pileup ..."

Rscript ${MODULE_DIR}/facets2n_snp-pileup-wrapper.R \
    --output-prefix "${FACETS2N_REFERENCE_DIR}/reference_normals" \
    --vcf-file "${DBSNP}" \
    --unmatched-normal-BAMS "${BAM_DIR}/*/*-N_recalibrated.bam"

## "==========================================================================="
## Reference loess normalization
## "==========================================================================="
echo "$(date +"%F") $(date +"%T") - Reference loess normalization ..."

Rscript - << 'EOF'
library(facets2n)

## Input and output files
pileup_file <- Sys.getenv("FACETS2N_REFERENCE_DIR")
pileup_file <- file.path(pileup_file, "reference_normals.snp_pileup.gz")
loess_file <- file.path(dirname(pileup_file), "reference_normals.loess.txt")

## Check if pileup file exists
if (!file.exists(pileup_file)) {
    stop(paste("Error: Pileup file not found:", pileup_file))
}

## Preprocess the snp-pileup file
pileup <- PreProcessSnpPileup(
    filename = pileup_file,
    is.Reference = TRUE
)

## Generate loess object
MakeLoessObject(
    pileup = pileup,
    write.loess = TRUE,
    outfilepath = loess_file,
    is_Reference = TRUE
)
EOF

## "==========================================================================="
## Generate the counts file
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
        out_dir="${FACETS2N_DIR}/${tumour_id}"
        rm -rf "${out_dir}"
        mkdir -p "${out_dir}"

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) - Generating counts file ..."
        ## Run FACETS snp-pileup
        Rscript ${MODULE_DIR}/facets2n_snp-pileup-wrapper.R \
            --vcf-file "${DBSNP}" \
            --normal-bam "${normal_bam}" \
            --tumor-bam "${tumour_bam}" \
            --unma
            --output-prefix "${out_dir}/${tumour_id}" \
            >& "${out_dir}/${tumour_id}.snp_pileup.log"

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