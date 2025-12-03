#!/bin/bash

#add GT to strelka vcf for downstream annovar
VAROUT=~/Strelka-Call/APLY/${TUMOUR}
echo $(date +"%F") $(date +"%T") "###########add GT to strelka vcf..................";
python ~/Scripts/DNA_analysis/add_GT_to_strelkaVCF.py $VAROUT/results/variants/${sample}_Strelka_normalized.vcf.gz $VAROUT/results/variants/${sample}_Strelka_normalized_GT.vcf
bgzip $VAROUT/results/variants/${sample}_Strelka_normalized_GT.vcf
tabix -p vcf $VAROUT/results/variants/${sample}_Strelka_normalized_GT.vcf.gz

  #Check header
  #bcftools query -l ./isec/0000.annot.vcf.gz
  #change strelka header
  bcftools reheader -s <(echo -e "SARC-006-N\nSARC-006-T") -o SARC-006_Strelka_reheader.vcf.gz SARC-006.strelka.ADVAF.vcf.gz
tabix SARC-006_Strelka_reheader.vcf.gz

# annotate AD (variant coverage)/VAF to strelka
java -Xmx4G -cp ~/Software/purple-2.26.jar com.hartwig.hmftools.purple.tools.AnnotateStrelkaWithAllelicDepth \
-in $VAROUT/results/variants/${sample}_Strelka_normalized_GT_reheader.vcf.gz \
-out $VAROUT/results/variants/${sample}.strelka.AD.vcf;

bcftools +fill-tags $VAROUT/results/variants/${sample}.strelka.AD.vcf -Ov -o $VAROUT/results/variants/${sample}.strelka.ADVAF.vcf -- -t FORMAT/VAF
bgzip $VAROUT/results/variants/${sample}.strelka.ADVAF.vcf
tabix $VAROUT/results/variants/${sample}.strelka.ADVAF.vcf.gz
rm $VAROUT/results/variants/${sample}.strelka.AD.vcf

# annotate AD (variant coverage) to mutect
bcftools +fill-tags SARC-004-T_normalized_filtered.vcf.gz -Ov -o SARC-006-T.mutect.AD.vcf -- -t FORMAT/VAF
bgzip SARC-006-T.mutect.AD.vcf
tabix SARC-006-T.mutect.AD.vcf.gz

  echo -e "##INFO=<ID=MuTect2,Number=0,Type=Flag,Description=\"Variant was called by MuTect2\">" > vcf.header
  echo -e "##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description=\"Variant failed filters in Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > vcf.rm.header
  echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" > vcf.map.header

# Get set differences of variant calls:
  # 0000: MuTect2 only
  # 0001: Strelka2 only
  # 0002: MuTect2 calls shared by Strelka2
  # 0003: Strelka2 calls shared by MuTect2

  bcftools isec \
    --output-type z \
    --prefix ./isec \
    /media/ezsharon/SF12T/DNA/SNV-call/Mutect/APYL/SARC-004-T/SARC-004-T.mutect.AD.vcf.gz /media/ezsharon/SF12T/DNA/SNV-call/Strelka/APYL/SARC-004-T/results/variants/SARC-004_Strelka_reheader.vcf.gz

bcftools annotate \
    --annotations ./isec/0003.vcf.gz \
    --include 'FILTER!="PASS"' \
    --mark-sites "+Strelka2FILTER" \
    -k \
    --output-type z \
    --output ./isec/0003.annot.vcf.gz \
    ./isec/0003.vcf.gz

  bcftools annotate \
    --header-lines /media/ezsharon/SF12T/DNA/SNV-call/merge/header/vcf.header \
    --annotations ./isec/0000.vcf.gz \
    --mark-sites +MuTect2 \
    --output-type z \
    --output ./isec/0000.annot.vcf.gz \
    ./isec/0000.vcf.gz

  bcftools annotate \
    --header-lines /media/ezsharon/SF12T/DNA/SNV-call/merge/header/vcf.header \
    --annotations ./isec/0002.vcf.gz \
    --mark-sites "+MuTect2;Strelka2" \
    --output-type z \
    --output ./isec/0002.tmp.vcf.gz \
    ./isec/0002.vcf.gz

  tabix --preset vcf ./isec/0002.tmp.vcf.gz
  tabix --preset vcf ./isec/0003.annot.vcf.gz

  bcftools annotate \
    --annotations ./isec/0003.annot.vcf.gz \
    --columns +INFO,+FORMAT,Strelka2FILTER \
    --output-type z \
    --output ./isec/0002.annot.vcf.gz \
    ./isec/0002.tmp.vcf.gz

  bcftools annotate \
    --header-lines /media/ezsharon/SF12T/DNA/SNV-call/merge/header/vcf.header \
    --annotations ./isec/0001.vcf.gz \
    --mark-sites +Strelka2 \
    --output-type z \
    --output ./isec/0001.annot.vcf.gz \
    ./isec/0001.vcf.gz

  tabix --preset vcf ./isec/0000.annot.vcf.gz
  tabix --preset vcf ./isec/0001.annot.vcf.gz
  tabix --preset vcf ./isec/0002.annot.vcf.gz

  # Concatenate the different sets
  bcftools concat \
    --allow-overlaps \
    --rm-dups all \
    ./isec/0000.annot.vcf.gz \
    ./isec/0001.annot.vcf.gz \
    ./isec/0002.annot.vcf.gz | \
  bcftools sort -Oz -o ./SARC-004-T.concatenated.sorted.vcf.gz
  
  #first filter
bcftools view -f PASS -Ov -o ./SARC-004-T.concatenated.pass.vcf.gz ./SARC-004-T.concatenated.sorted.vcf.gz
bcftools view -i 'FILTER!="PASS"' -Ov -o ./SARC-004-T.concatenated.unpass.vcf.gz ./SARC-004-T.concatenated.sorted.vcf.gz

#annotate with repeatmasker and blacklists
bcftools annotate \
  --header-lines /media/ezsharon/SF12T/DNA/SNV-call/merge/header/vcf.rm.header \
  --annotations ~/Reference/RepeatMasker.bed.gz \
  --columns CHROM,FROM,TO,RepeatMasker \
  ./SARC-004-T.concatenated.pass.vcf.gz | \
bcftools annotate \
  --header-lines /media/ezsharon/SF12T/DNA/SNV-call/merge/header/vcf.map.header \
  --annotations ~/Reference/blacklist.bed.gz \
  --columns CHROM,FROM,TO,EncodeDacMapability \
  /dev/stdin \
  --output-type z \
  --output ./SARC-004-T.union.vcf.gz
    #repeatmasker: https://genome.ucsc.edu/cgi-bin/hgTables. variant falls within a region of the genome that is known to contain repetitive sequences. Variants in these regions are often treated with caution because the repetitive nature of the DNA can lead to alignment errors and false variant calls.
    #blacklist: https://www.encodeproject.org/files/ENCFF269URO/. genomic regions that are known to be problematic for variant calling, often filtered out or ignored in downstream analyses.

tabix --preset vcf ./SARC-004-T.union.vcf.gz

# Filter out variants in RepeatMasker or Mapability
bcftools filter -e 'INFO/RepeatMasker != "." || INFO/EncodeDacMapability != "."' -Oz -o ./SARC-004-T.union.filtered.vcf.gz ./SARC-004-T.union.vcf.gz


#Annotate by Annovar
echo $(date +"%F") $(date +"%T") "###########Annotating by annovar..................";
perl ~/Software/annovar/table_annovar.pl ./SARC-004-T.union.vcf.gz \
~/Software/annovar/humandb/ \
-buildver hg38 -out ./SARC-004-T.union \
-remove \
-protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
-operation gx,r,f,f,f,f,f \
-nastring . -polish -xreffile ~/Software/annovar/example/gene_fullxref.txt \
--otherinfo --vcfinput ;


less -S ./SARC-004-T.union.hg38_multianno.txt | awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "RAD_N", "VAD_N", "DP_N", "VAF_N", "RAD_T", "VAD_T", "DP_T", "VAF_T"} NR >1 {split($(NF-1), a, ":"); split(a[2], ad, ","); split($(NF), b, ":"); split(b[2], bd, ",");  $(NF+1)=ad[1]; $(NF+2)=ad[2]; $(NF+3)=ad[1]+ad[2]; if (length(a) == 11) $(NF+4)=a[11]; else if (length(a) == 16) $(NF+4)=a[9]; $(NF+5)=bd[1]; $(NF+6)=bd[2]; $(NF+7)=bd[1]+bd[2]; if (length(b) == 11) $(NF+8)=b[11]; else if (length(b) == 16) $(NF+8)=b[9]; print $0} ' > ./SARC-004-T.union.txt


