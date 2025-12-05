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
export REFERENCE_DIR="/lustre1/g/path_my/Reference"
export SOFTWARE_DIR="/lustre1/g/path_my/Software"

export CONTAINER_DIR="${PROJECT_DIR}/containers"
export MODULES_DIR="${PROJECT_DIR}/scripts/modules"
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed.gz"
export BAM_DIR="${PROJECT_DIR}/data/SARC/BAM"
export SV_CALL_DIR="${PROJECT_DIR}/data/SARC/SV-call"
export STRELKA_DIR="${PROJECT_DIR}/data/SARC/Strelka2"
mkdir -p "${STRELKA_DIR}"

PCGR_REFERENCE="${REFERENCE_DIR}/PCGR_reference/20250314"
VEP_CACHE="${REFERENCE_DIR}/VEP_cache"

# export TUMOR_DP_MIN=20
# export TUMOR_AF_MIN=0.05
# export CONTROL_DP_MIN=10
# export CONTROL_AF_MAX=0.01
# export TMB_DP_MIN=20
# export TMB_AF_MIN=0.05
# export N_COPY_GAIN=3

## "==========================================================================="
## Function to run strelka for a sample
## "==========================================================================="
strelka() {
    
    local tumour_id="$1"
    
    tumour_id="SARC-035-T"
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

        ## "==================================================================="
        ## Run Strelka configuration
        ## "==================================================================="
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

        ## "==================================================================="
        ## Run the Strelka workflow
        ## "==================================================================="
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Strelka workflow ..."
        "${out_dir}/runWorkflow.py" \
            -m local \
            -j 16 \
            >& "${out_dir}/strelka_workflow.log"

        ## "==================================================================="
        ## Combine the SNV and Indel VCFs
        ## "==================================================================="
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

        ## "==================================================================="
        ## Annovar annotation
        ## "==================================================================="
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Annovar annotation ..."
        ## Add GT to Strelka VCF for downstream annovar
        singularity exec \
            --bind "${PROJECT_DIR}":"${PROJECT_DIR}" \
            "${CONTAINER_DIR}/pysam.sif" \
            python "${MODULES_DIR}/add_GT_to_strelkaVCF.py" \
                "${out_dir}/results/variants/${tumour_id}_strelka_normalized.vcf.gz" \
                "${out_dir}/results/variants/${tumour_id}_strelka_normalized_GT.vcf"
        
        bgzip "${out_dir}/results/variants/${tumour_id}_strelka_normalized_GT.vcf"
        tabix -p vcf "${out_dir}/results/variants/${tumour_id}_strelka_normalized_GT.vcf.gz"

        annovar_output_dir="${out_dir}/Annovar"
        mkdir -p "${annovar_output_dir}"
        
        ## Run Annovar annotation
        perl "${SOFTWARE_DIR}/annovar/table_annovar.pl" \
            "${out_dir}/results/variants/${tumour_id}_strelka_normalized_GT.vcf.gz" \
            "${SOFTWARE_DIR}/annovar/humandb/" \
            -buildver hg38 \
            -out "${annovar_output_dir}/${tumour_id}" \
            -remove \
            -protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
            -operation gx,r,f,f,f,f,f \
            -nastring . \
            -polish \
            -xreffile "${REFERENCE_DIR}/annovar/example/gene_fullxref.txt" \
            --otherinfo \
            --vcfinput \
            >& "${annovar_output_dir}/${tumour_id}.annovar.log"

        less -S "${annovar_output_dir}/${tumour_id}.hg38_multianno.txt" |
            awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "AD", "AF", "DP"}; NR >1 {split($NF, a, ":"); $(NF+1)=a[2]; $(NF+1)=a[3]; $(NF+1)=a[4]; print}' > \
                "${annovar_output_dir}/${tumour_id}.annovar.txt"
        
        ## "==================================================================="
        ## PCGR annotation
        ## "==================================================================="
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) PCGR annotation ..."
        pcgr_output_dir="${out_dir}/PCGR"
        rm -rf "${pcgr_output_dir}"
        mkdir -p "${pcgr_output_dir}"

        ## Reformat VCF for PCGR
        reformatted_vcf="${pcgr_output_dir}/${tumour_id}.reformatted.vcf.gz"
        singularity exec \
            --bind "${out_dir}/results/variants":"${out_dir}/results/variants" \
            --bind "${pcgr_output_dir}":"${pcgr_output_dir}" \
            --bind "${MODULES_DIR}":"${MODULES_DIR}" \
            --bind "/tmp":"/tmp" \
            "${CONTAINER_DIR}/pysam.sif" \
            python "${MODULES_DIR}/pcgr_reformat_vcf_strelka_tumour_normal.py" \
                -I "${out_dir}/results/variants/${tumour_id}_strelka_normalized_GT.vcf.gz" \
                -O "${reformatted_vcf}"
        
        tabix -p vcf "${reformatted_vcf}"

        ## Run PCGR annotation
        singularity exec \
            --bind "${PROJECT_DIR}":"${PROJECT_DIR}" \
            --bind "${PCGR_REFERENCE}":"${PCGR_REFERENCE}" \
            --bind "${VEP_CACHE}":"${VEP_CACHE}" \
            --bind "${pcgr_output_dir}":"${pcgr_output_dir}" \
            --bind "${out_dir}/results/variants":"${out_dir}/results/variants" \
            "${CONTAINER_DIR}/pcgr.sif" \
                pcgr \
                --input_vcf "${reformatted_vcf}" \
                --vep_dir "${VEP_CACHE}" \
                --refdata_dir "${PCGR_REFERENCE}" \
                --output_dir "${pcgr_output_dir}" \
                --genome_assembly grch38 \
                --sample_id "${tumour_id}" \
                --assay WES \
                --effective_target_size_mb 34 \
                --tumor_dp_tag TDP \
                --tumor_af_tag TAF \
                --control_dp_tag NDP \
                --control_af_tag NAF \
                --tumor_dp_min "${TUMOR_DP_MIN}" \
                --tumor_af_min "${TUMOR_AF_MIN}" \
                --estimate_tmb \
                --tmb_dp_min "${TMB_DP_MIN}" \
                --tmb_af_min "${TMB_AF_MIN}" \
                --estimate_msi \
                --estimate_signatures \
                --vcf2maf \
                --ignore_noncoding \
                --force_overwrite \
                >& "${pcgr_output_dir}/${tumour_id}.pcgr.log"

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