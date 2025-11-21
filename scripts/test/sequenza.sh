#!/bin/bash

# Test the singularity container for Sequenza
singularity shell \
    --bind /home/zhonggr/projects/250224_sarcoma_multiomics/test:/data \
    --bind /home/zhonggr/projects/250224_sarcoma_multiomics/data/reference:/reference \
    /home/zhonggr/projects/250224_sarcoma_multiomics/containers/singularity/sequenzautils.sif

export REFERENCE=/reference/Gencode/gencode.hg38.v36.primary_assembly.fa
export GC=/reference/Gencode/hg38.gc50Base.wig.gz
export tumour=/data/DFSP-001-T/DFSP-001-T_recalibrated.bam
export normal=/data/DFSP-001-N/DFSP-001-N_recalibrated.bam
sequenza_dir=/data/Sequenza

# export REFERENCE=/home/zhonggr/projects/250224_sarcoma_multiomics/data/reference/Gencode/gencode.hg38.v36.primary_assembly.fa
# export GC=/home/zhonggr/projects/250224_sarcoma_multiomics/data/reference/Gencode/hg38.gc50Base.wig.gz
# export tumour=/home/zhonggr/projects/250224_sarcoma_multiomics/test/DFSP-001-T/DFSP-001-T_recalibrated.bam
# export normal=/home/zhonggr/projects/250224_sarcoma_multiomics/test/DFSP-001-N/DFSP-001-N_recalibrated.bam
# export sequenza_dir=/home/zhonggr/projects/250224_sarcoma_multiomics/test/Sequenza

mkdir -p ${sequenza_dir}

# Check if reference, BAM files, and GC file exist
if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference file ${REFERENCE} not found"
    exit 1
fi

if [ ! -f "${GC}" ]; then
    echo "ERROR: GC content file ${GC} not found"
    exit 1
fi

if [ ! -f "${tumour}" ]; then
    echo "ERROR: Tumor BAM file ${tumour} not found"
    exit 1
fi

if [ ! -f "${normal}" ]; then
    echo "ERROR: Normal BAM file ${normal} not found"
    exit 1
fi

sequenza-utils bam2seqz \
    --normal ${normal} \
    --tumor ${tumour} \
    -gc ${GC} \
    --fasta ${REFERENCE} \
    --output ${sequenza_dir}/DFSP-001-T.seqz.gz

sequenza-utils seqz_binning \
    -s ${sequenza_dir}/DFSP-001-T.seqz.gz \
    -w 50 \
    -o ${sequenza_dir}/DFSP-001-T.bin50.seqz.gz



