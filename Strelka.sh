#!/bin/bash

SERVER=/run/user/1000/gvfs/smb-share:server=maximus-nas.local,share=genomics
REFERENCE=/media/Data/Reference/Gencode/gencode.hg38.v36.primary_assembly.fa
INTERVAL=/media/Data/Reference/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed.gz
BAM=/media/Data/DFSP/Test

CASELIST=$(ls $BAM | cut -d '-' -f1,2 | uniq)

for sample in $CASELIST; do 
echo $sample

if [ $sample == "DFSP-PF1" ] || [ $sample == "DFSP-PF2" ] || [ $sample == "DFSP-PF3" ];
then
TUMOUR=${sample}-T;
NORMAL=$(echo $sample | sed 's/PF/PC/')-N;
else
TUMOUR=${sample}-T;
NORMAL=${sample}-N;

fi

VAROUT=/media/Data/DFSP/Strelka-Call/${sample}
echo $TUMOUR
echo $NORMAL

/media/Data/Software/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
--normalBam $BAM/$NORMAL/${NORMAL}_recalibrated.bam \
--tumorBam $BAM/$TUMOUR/${TUMOUR}_recalibrated.bam \
--callRegions $INTERVAL \
--exome \
--referenceFasta $REFERENCE \
--indelCandidates /media/Data/DFSP/Manta-call/${sample}/results/variants/candidateSmallIndels.vcf.gz \
--runDir $VAROUT

$VAROUT/runWorkflow.py -m local -j 16

bcftools concat -Oz -a $VAROUT/results/variants/somatic.snvs.vcf.gz $VAROUT/results/variants/somatic.indels.vcf.gz -o $VAROUT/results/variants/${sample}_Strelka.vcf.gz
bcftools norm -m-both -f $REFERENCE -Oz -o $VAROUT/results/variants/${sample}_Strelka_normalized.vcf.gz $VAROUT/results/variants/${sample}_Strelka.vcf.gz
tabix $VAROUT/results/variants/${sample}_Strelka_normalized.vcf.gz

done
