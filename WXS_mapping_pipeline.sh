#!/bin/bash

SAMPLENAME=$1

WXDIR=/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/AML/exomes/cell_lines/fastq/

BWAINDEX=/srv/shared/vanloo/pipeline-files/human/references/alignment/GRCh37-1000G/bwa_index/genome.fa
REFGENOME=/srv/shared/vanloo/pipeline-files/human/references/alignment/GRCh37-1000G/genome.fa
TARGETPANEL=/srv/shared/vanloo/pipeline-files/human/references/targeted_panels/SureSelect_Human_All_Exon_V6_UTR_r2/S07604624_Padded_GRCh37.bed
KNOWNINDELS=/srv/shared/vanloo/pipeline-files/human/references/b37indels_and_snps/Mills_and_1000G_gold_standard.indels.b37.vcf.bgz


echo "Running on sample ${SAMPLENAME}"

FWDREADS=$(find ${WXDIR} -name "${SAMPLENAME}*_R1_001.fastq.gz")
REVREADS=$(find ${WXDIR} -name "${SAMPLENAME}*_R3_001.fastq.gz")

CWD=$(pwd)

OUTDIR=/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/AML/exomes/cell_lines/mapped/${SAMPLENAME}_mapped

mkdir -p ${OUTDIR}
cd ${OUTDIR}


echo "Running trimmomatic"
# module purge
# ml Trimmomatic

java -Xmx48G -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.38.jar \
      PE -threads 8 -phred33 \
      ${FWDREADS} ${REVREADS} \
      ${SAMPLENAME}_R1_trimmed.fastq.gz ${SAMPLENAME}_R1_trimmed_unpaired.fastq.gz \
      ${SAMPLENAME}_R3_trimmed.fastq.gz ${SAMPLENAME}_R3_trimmed_unpaired.fastq.gz \
      ILLUMINACLIP:${EBROOTTRIMMOMATIC}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


echo "Running bwa-mem"
# ml BWA
bwa mem -M -a -t 8 -R "@RG\tID:${SAMPLENAME}\tPL:ILLUMINA\tLB:${SAMPLENAME}\tSM:${SAMPLENAME}" \
      ${BWAINDEX} \
      ${SAMPLENAME}_R1_trimmed.fastq.gz ${SAMPLENAME}_R3_trimmed.fastq.gz \
      > ${SAMPLENAME}.sam


echo "Running samtools/picard"
# module purge
# ml SAMtools
# ml picard

samtools sort -@ 8 -m 6G -o ${SAMPLENAME}_sorted.bam ${SAMPLENAME}.sam

samtools index ${SAMPLENAME}_sorted.bam

java -Xmx48G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=${SAMPLENAME}_sorted.bam \
      O=${SAMPLENAME}_sorted_unique.bam \
      M=${SAMPLENAME}_marked_dup_metrics.txt \
      ASSUME_SORTED=true

samtools index ${SAMPLENAME}_sorted_unique.bam


echo "Running GATK indel realignment"
# module purge
# ml GATK/3.8-1-Java-1.8.0_162

java -Xmx48G -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator \
      -nt 6 \
      -R ${REFGENOME} \
      -L ${TARGETPANEL} \
      --known ${KNOWNINDELS} \
      -I ${SAMPLENAME}_sorted_unique.bam \
      -o ${SAMPLENAME}_sorted_unique.bam.list

java -Xmx48G -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner \
      -targetIntervals ${SAMPLENAME}_sorted_unique.bam.list \
      -R ${REFGENOME} \
      --knownAlleles ${KNOWNINDELS} \
      -I ${SAMPLENAME}_sorted_unique.bam \
      -o ${SAMPLENAME}_sorted_unique_realigned.bam \
      --filter_bases_not_stored


echo "Cleaning up"

# module purge
rm ${SAMPLENAME}*_trimmed.fastq.gz ${SAMPLENAME}*_trimmed_unpaired.fastq.gz \
    ${SAMPLENAME}.sam ${SAMPLENAME}.bam ${SAMPLENAME}_sorted.bam* \
    ${SAMPLENAME}_sorted_unique.bam* ${SAMPLENAME}_marked_dup_metrics.txt

cd ${CWD}
