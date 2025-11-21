#!/bin/bash
#SBATCH --job-name=SomaticSV_call
#SBATCH --partition=amd
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.out
#SBATCH --error=/lustre1/g/path_my/pipeline/somatic_variants_calling/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk


# Parallel Somatic SV calling using Delly + Manta + AnnotSV
# This refactor runs one tumour_id sample per job using GNU Parallel.

export PROJECT_DIR=/lustre1/g/path_my/pipeline/somatic_variants_calling
export REFERENCE_DIR=/lustre1/g/path_my/Reference
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export REGIONS="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed.gz"

export SAMTOOLS_PATH=/home/zhonggr/miniforge3/bin/samtools
export CONTAINER_DIR=/lustre1/g/path_my/pipeline/somatic_variants_calling/containers

export BAM_DIR=${PROJECT_DIR}/data/DFSP/BAM

export SV_CALL_DIR=${PROJECT_DIR}/data/DFSP/SV-call

export SOFTWARE_DIR=/lustre1/g/path_my/Software
export EXCLUDE="${SOFTWARE_DIR}/delly/excludeTemplates/human.hg38.excl.tsv"
export ANNOTSV_PATH=${SOFTWARE_DIR}/AnnotSV

run_somatic_sv() {
    
    tumour_id="$1"

    patient_id=$(echo "$tumour_id" | cut -d'-' -f1,2)
    normal_id="${patient_id}-N"

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

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Delly Calling ..."
        
        delly_dir="${SV_CALL_DIR}/${tumour_id}/Delly"
        rm -rf "${delly_dir}"
        mkdir -p "${delly_dir}"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${SOFTWARE_DIR}:${SOFTWARE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/delly.sif" \
            delly call \
                -x "${EXCLUDE}" \
                -o "${delly_dir}/${tumour_id}.bcf" \
                -g "$REFERENCE" \
                "${tumour_bam}" \
                "${normal_bam}" \
                >& "${delly_dir}/${tumour_id}.delly.call.log"

        echo -e "${tumour_id}\ttumor\n${normal_id}\tcontrol" > "${delly_dir}/samples.tsv"

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Delly Filtering ..."
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/delly.sif" \
            delly filter \
                -f somatic \
                -o "${delly_dir}/${tumour_id}.somatic.filter.bcf" \
                -s "${delly_dir}/samples.tsv" \
                "${delly_dir}/${tumour_id}.bcf"
        
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/bcftools.sif" \
            bcftools view "${delly_dir}/${tumour_id}.somatic.filter.bcf" -Oz > "${delly_dir}/${tumour_id}.somaticSV.delly.vcf.gz"
        
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/tabix.sif" \
            tabix "${delly_dir}/${tumour_id}.somaticSV.delly.vcf.gz"

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Manta Config ..."
        manta_dir="${SV_CALL_DIR}/${tumour_id}/Manta"
        rm -rf "${manta_dir}"
        mkdir -p "${manta_dir}"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/manta.sif" \
            configManta.py \
                --normalBam="${normal_bam}" \
                --tumourBam="${tumour_bam}" \
                --exome \
                --referenceFasta="${REFERENCE}" \
                --runDir="${manta_dir}" \
                --callRegions="${REGIONS}" \
                >& "${manta_dir}/${tumour_id}.manta.config.log"
        
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Manta Calling ..."
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/manta.sif" \
            "$manta_dir/runWorkflow.py" -j 8 \
            >& "${manta_dir}/${tumour_id}.manta.runWorkflow.log"
        
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Manta Conversion ..."
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind "${SAMTOOLS_PATH}:${SAMTOOLS_PATH}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/manta.sif" \
            python2 /opt/conda/libexec/convertInversion.py \
                "${SAMTOOLS_PATH}" \
                "$REFERENCE" \
                "$manta_dir/results/variants/somaticSV.vcf.gz" \
                >& "${manta_dir}/${tumour_id}.manta.convert.log"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/tabix.sif" \
            tabix -f --preset vcf "$manta_dir/results/variants/somaticSV.vcf.gz"

        ## Reorganize outputs
        mv ${manta_dir}/results/variants/candidateSmallIndels.vcf.gz \
        ${manta_dir}/${tumour_id}.candidateSmallIndels.vcf.gz
        mv ${manta_dir}/results/variants/candidateSmallIndels.vcf.gz.tbi \
        ${manta_dir}/${tumour_id}.candidateSmallIndels.vcf.gz.tbi
        mv ${manta_dir}/results/variants/candidateSV.vcf.gz \
        ${manta_dir}/${tumour_id}.candidateSV.vcf.gz
        mv ${manta_dir}/results/variants/candidateSV.vcf.gz.tbi \
        ${manta_dir}/${tumour_id}.candidateSV.vcf.gz.tbi
        mv ${manta_dir}/results/variants/diploidSV.vcf.gz \
        ${manta_dir}/${tumour_id}.diploidSV.vcf.gz
        mv ${manta_dir}/results/variants/diploidSV.vcf.gz.tbi \
        ${manta_dir}/${tumour_id}.diploidSV.vcf.gz.tbi
        mv ${manta_dir}/results/variants/somaticSV.vcf.gz \
        ${manta_dir}/${tumour_id}.somaticSV.manta.vcf.gz
        mv ${manta_dir}/results/variants/somaticSV.vcf.gz.tbi \
        ${manta_dir}/${tumour_id}.somaticSV.manta.vcf.gz.tbi

        ## Merge outputs from Delly and Manta
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Merging Delly and Manta outputs ..."
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/bcftools.sif" \
            bcftools view \
            --samples "${tumour_id},${normal_id}" \
            --output-type z \
            --output-file "${manta_dir}/${tumour_id}.manta.swap.vcf.gz" \
            "${manta_dir}/${tumour_id}.somaticSV.manta.vcf.gz"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/tabix.sif" \
            tabix --preset vcf "${manta_dir}/${tumour_id}.manta.swap.vcf.gz"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/bcftools.sif" \
            bcftools concat \
            --allow-overlaps \
            --output-type z \
            --output "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.delly.manta.unfiltered.vcf.gz" \
            "${delly_dir}/${tumour_id}.somaticSV.delly.vcf.gz" \
            "${manta_dir}/${tumour_id}.manta.swap.vcf.gz" \
            >& "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.bcftools.concat.log"
        
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/tabix.sif" \
            tabix --preset vcf "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.delly.manta.unfiltered.vcf.gz"

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Filtering merged SV ..."
        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/bcftools.sif" \
            bcftools filter \
            --include 'FILTER="PASS"' \
            --output-type z \
            --output "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.delly.manta.filtered.vcf.gz" \
            "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.delly.manta.unfiltered.vcf.gz"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/bcftools.sif" \
            bcftools sort \
            --output-type z \
            --output-file "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.delly.manta.vcf.gz" \
            "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.delly.manta.filtered.vcf.gz"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${SV_CALL_DIR}:${SV_CALL_DIR}" \
            --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
            --bind "${BAM_DIR}:${BAM_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/tabix.sif" \
            tabix --preset vcf "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.delly.manta.vcf.gz"

        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) AnnotSV ..."
        ${ANNOTSV_PATH}/bin/AnnotSV \
            -SVinputFile "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.delly.manta.vcf.gz" \
            -outputFile "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.somaticSV.annotated.tsv" \
            -genomeBuild GRCh38 \
            -annotationMode both \
            -SVminSize 50 \
            >& "${SV_CALL_DIR}/${tumour_id}/${tumour_id}.AnnotSV.log"
    fi

}

export -f run_somatic_sv

# Collect tumour_id ids (directories in BAM that contain 'T' in the name)
# tumour_ids=$(find "${BAM_DIR}" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | grep "T" | sort)
# PARALLEL_JOBS=20

# echo "${tumour_ids}" | parallel --jobs "${PARALLEL_JOBS}" run_somatic_sv {}


## Find out samples that do not have correct results files in SV_CALL_DIR

# Get all tumour samples from BAM directory (directories containing 'T' in the name)
tumour_ids=$(find "${BAM_DIR}" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | grep "T" | sort)
# tumour_ids=$(find "${SV_CALL_DIR}" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | grep "T" | sort)

## Find samples that have completed SV calling results
samples_with_results=""
for sample in $tumour_ids; do
    if [[ -f "${SV_CALL_DIR}/${sample}/${sample}.somaticSV.annotated.tsv" ]]; then
        samples_with_results="${samples_with_results}${sample}\n"
    fi
done

## Number of samples with results
num_samples_with_results=$(echo -e "$samples_with_results" | grep -v '^$' | wc -l)

# Get samples that don't have SV results files but have the required BAM files
samples_missing_results=""
for sample in $tumour_ids; do
    patient_id=$(echo "$sample" | cut -d'-' -f1,2)
    normal_id="${patient_id}-N"
    
    # Check if both tumour and normal BAM files exist
    tumour_bam="${BAM_DIR}/${sample}/${sample}_recalibrated.bam"
    normal_bam="${BAM_DIR}/${normal_id}/${normal_id}_recalibrated.bam"
    
    if [[ -f "${tumour_bam}" ]] && [[ -f "${normal_bam}" ]] && [[ ! -f "${SV_CALL_DIR}/${sample}/${sample}.somaticSV.annotated.tsv" ]]; then
        samples_missing_results="${samples_missing_results}${sample}\n"
    fi
done
samples_missing_results=$(echo -e "$samples_missing_results" | grep -v '^$')

echo "Found $(echo -e "$samples_missing_results" | wc -l) samples missing SV calling results files:"
echo -e "$samples_missing_results"

# Run SV calling for missing samples
PARALLEL_JOBS=2
if [[ -n "$samples_missing_results" ]]; then
    echo -e "$samples_missing_results" | parallel --jobs "$PARALLEL_JOBS" run_somatic_sv {}
fi
