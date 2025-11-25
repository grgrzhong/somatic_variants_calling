#!/bin/bash
#SBATCH --job-name=cnvkit
#SBATCH --partition=amd
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_Multiomics/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_Multiomics/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

##############################################################################
# Script to run CNVkit for CNV analysis on tumor samples
# https://cnvkit.readthedocs.io/en/stable/quickstart.html
##############################################################################
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cnvkit

## Define directories and reference files
ref_dir=/home/zhonggr/projects/250224_sarcoma_multiomics/data/reference
reference=${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa
target=${ref_dir}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed
bam_dir=/home/zhonggr/projects/250224_sarcoma_multiomics/data/WES/bam

## Define output directory
work_dir=/home/zhonggr/projects/250224_sarcoma_multiomics/data/WES/cnvkit
mkdir -p ${work_dir}
cnvkit_ref_dir=${work_dir}/cnvkit_reference
cnvkit_norm_ref=${cnvkit_ref_dir}/cnvkit_pooled_normal_reference.cnn

## Download the gene annotations file as the targets not provided
# wget -P ${cnvkit_ref_dir} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
# gunzip -f ${cnvkit_ref_dir}/refFlat.txt.gz

## Collect the normal bam files
normal_bams=$(find ${bam_dir} -type f -name "*.bam" | grep "N")

## Step1: Create a pooled reference of per-bin copy number estimates from several normal samples
if [ ! -f "${cnvkit_norm_ref}" ]; then

    echo "Creating pooled reference..."

    cnvkit.py batch \
        --normal ${normal_bams} \
        --fasta ${reference} \
        --targets ${target} \
        --processes 10 \
        --annotate ${cnvkit_ref_dir}/refFlat.txt \
        --output-reference ${cnvkit_norm_ref} \
        --output-dir ${cnvkit_ref_dir} \
        >& ${work_dir}/cnvkit_pooled_normal_reference.log
fi

## Step2. Use this reference in processing all tumor samples
tumour_bams=$(find ${bam_dir} -type f -name "*.bam" | grep -v "N")

# run_cnvkit() {
#     local tumour_bam=$1
#     local work_dir=$2
#     local cnvkit_norm_ref=$3
    
#     tumour_id=$(basename ${tumour_bam} _recalibrated.bam) 
#     out_dir=${work_dir}/${tumour_id}

#     echo "Processing : ${tumour_id}"

#     # Run cnvkit for the tumour sample
#     cnvkit.py batch \
#         ${tumour_bam} \
#         --reference ${cnvkit_norm_ref} \
#         --scatter \
#         --diagram \
#         --output-dir ${out_dir} \
#         >& ${out_dir}/${tumour_id}.cnvkit.log
    
#     # Generate additional reports
#     # cnvkit.py heatmap \
#     #     "${out_dir}/${tumour_id}.cns" \
#     #     -o "${sample_output_dir}/${sample_id}_heatmap.pdf"
        
#     # # Generate genome-wide plot
#     # cnvkit.py scatter "${out_dir}/${tumour_id}.cnr" \
#     #     -s "${out_dir}/${tumour_id}.cns" \
#     #     -o "${out_dir}/${tumour_id}.scatter.pdf"

# }

# echo "${tumor_bams}" | parallel \
#     -j 4 \
#     run_cnvkit {} "${work_dir}" "${cnvkit_norm_ref}"

# echo "All tumor samples processed successfully!"