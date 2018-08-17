#!/bin/bash

SAMPLENAME=$1
WXDIR=/camp/lab/vanloop/working/demeulj/projects/2016_mansour_ASE_T-ALL/data/Cell_lines/Exomes_TALL_Cell_Lines/

OUTDIR=${WXDIR}${SAMPLENAME}
FWDREADS=$(find ${WXDIR} -name "${SAMPLENAME}*_R1_001.fastq.gz")
REVREADS=$(find ${WXDIR} -name "${SAMPLENAME}*_R3_001.fastq.gz")

BWAREFGENOME=/camp/lab/vanloop/reference/Genomics/babs/homo_sapiens/ensembl/GRCh37/release-75/genome_idx/bwa/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa
PICARD=/home/camp/demeulj/.local/easybuild/software/picard/2.17.8-Java-1.8.0_131/picard.jar
GATK=/home/camp/demeulj/gatk-4.0.1.1/gatk


echo "Running on sample ${SAMPLENAME}"
mkdir ${OUTDIR}
cd ${OUTDIR}


echo "Running trimmomatic"
module purge
ml Trimmomatic

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 4 -phred33 \
${FWDREADS} ${REVREADS} \
${SAMPLENAME}_R1_trimmed.fastq.gz ${SAMPLENAME}_R1_trimmed_unpaired.fastq.gz \
${SAMPLENAME}_R2_trimmed.fastq.gz ${SAMPLENAME}_R2_trimmed_unpaired.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


echo "Running bwa-mem"
ml BWA
bwa mem -t 4 ${BWAREFGENOME} ${SAMPLENAME}_R1_trimmed.fastq.gz ${SAMPLENAME}_R2_trimmed.fastq.gz > ${SAMPLENAME}.sam


echo "Running samtools/picard"
module purge
ml SAMtools
ml Java

samtools view -@ 4 -b ${SAMPLENAME}.sam -o ${SAMPLENAME}.bam

java -Xmx24G -jar ${PICARD} SortSam \
      I=${SAMPLENAME}.bam \
      O=${SAMPLENAME}_qnamesort.bam \
      SORT_ORDER=queryname

java -Xmx24G -jar ${PICARD} MarkDuplicates \
      I=${SAMPLENAME}_qnamesort.bam \
      O=${SAMPLENAME}_qnamesort_markdup.bam \
      M=${SAMPLENAME}_marked_dup_metrics.txt \
      ASSUME_SORT_ORDER=queryname

java -Xmx24G -jar ${PICARD} AddOrReplaceReadGroups \
      I=${SAMPLENAME}_qnamesort_markdup.bam \
      O=${SAMPLENAME}_clean.bam \
      SORT_ORDER=coordinate \
      RGID=1 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=${SAMPLENAME}

samtools index ${SAMPLENAME}_clean.bam


echo "Cleaning up"

module purge
rm ${SAMPLENAME}*_trimmed.fastq.gz ${SAMPLENAME}*_trimmed_unpaired.fastq.gz \
    ${SAMPLENAME}.sam ${SAMPLENAME}.bam ${SAMPLENAME}_qnamesort.bam \
    ${SAMPLENAME}_qnamesort_markdup.bam ${SAMPLENAME}_marked_dup_metrics.txt

cd $WXDIR
