#!/bin/bash

SAMPLENAME=$1
RNADIR=/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/Cell_lines/RNA_Seq_TALL_Cell_Lines/

OUTDIR=${RNADIR}${SAMPLENAME}_mapped

STARREFGENOME=/srv/shared/vanloo/pipeline-files/human/references/alignment/hg19/STARidx/
STARREFANNOT=/srv/shared/vanloo/pipeline-files/human/references/annotation/GENCODE/gencode.v23lift37.annotation.gtf

mkdir ${OUTDIR}
cd ${OUTDIR}

module purge
ml Trimmomatic

# for lane in $(seq 1 4)
# do
#       FWDREADS=$(find ${RNADIR} -name "*${SAMPLENAME}*_L00${lane}_R1_001.fastq.gz")
#       REVREADS=$(find ${RNADIR} -name "*${SAMPLENAME}*_L00${lane}_R2_001.fastq.gz")

#       echo "Running trimmomatic on sample ${SAMPLENAME}, lane ${lane}"

#       java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 4 -phred33 \
#             ${FWDREADS} ${REVREADS} \
#             ${SAMPLENAME}_L00${lane}_R1_trimmed.fastq.gz ${SAMPLENAME}_L00${lane}_R1_trimmed_unpaired.fastq.gz \
#             ${SAMPLENAME}_L00${lane}_R2_trimmed.fastq.gz ${SAMPLENAME}_L00${lane}_R2_trimmed_unpaired.fastq.gz \
#             ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# done

trim () {
      local lane=$1
      FWDREADS=$(find ${RNADIR} -name "*${SAMPLENAME}*_L00${lane}_R1_001.fastq.gz")
      REVREADS=$(find ${RNADIR} -name "*${SAMPLENAME}*_L00${lane}_R2_001.fastq.gz")

      java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 4 -phred33 \
            ${FWDREADS} ${REVREADS} \
            ${SAMPLENAME}_L00${lane}_R1_trimmed.fastq.gz ${SAMPLENAME}_L00${lane}_R1_trimmed_unpaired.fastq.gz \
            ${SAMPLENAME}_L00${lane}_R2_trimmed.fastq.gz ${SAMPLENAME}_L00${lane}_R2_trimmed_unpaired.fastq.gz \
            ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
}

# for run in $(seq 1 4); do trim ${run} & done
# wait

# echo "Running STAR"
# ml BWA
# bwa mem -t 4 ${BWAREFGENOME} ${SAMPLENAME}_R1_trimmed.fastq.gz ${SAMPLENAME}_R2_trimmed.fastq.gz > ${SAMPLENAME}.sam

# module purge
# ml STAR


# STAR --runThreadN 12 \
#       --genomeDir ${STARREFGENOME} \
#       --readFilesCommand zcat \
#       --readFilesIn ${SAMPLENAME}_L001_R1_trimmed.fastq.gz,${SAMPLENAME}_L002_R1_trimmed.fastq.gz,${SAMPLENAME}_L003_R1_trimmed.fastq.gz,${SAMPLENAME}_L004_R1_trimmed.fastq.gz ${SAMPLENAME}_L001_R2_trimmed.fastq.gz,${SAMPLENAME}_L002_R2_trimmed.fastq.gz,${SAMPLENAME}_L003_R2_trimmed.fastq.gz,${SAMPLENAME}_L004_R2_trimmed.fastq.gz \
#       --outFileNamePrefix ${SAMPLENAME}_ \
#       --outSAMtype BAM Unsorted \
#       --outSAMattrRGline ID:L001 SM:${SAMPLENAME} LB:LIB1 , ID:L002 SM:${SAMPLENAME} LB:LIB1 , ID:L003 SM:${SAMPLENAME} LB:LIB1 , ID:L004 SM:${SAMPLENAME} LB:LIB1 \
#       --outSAMattributes NH NM MD \
#       --outSAMunmapped Within \
#       --outSAMmapqUnique 50 \
#       --outSJfilterCountUniqueMin -1 2 2 2 \
#       --outSJfilterCountTotalMin -1 2 2 2 \
#       --outFilterType BySJout \
#       --outFilterIntronMotifs RemoveNoncanonical \
#       --chimSegmentMin 12 \
#       --chimScoreDropMax 30 \
#       --chimScoreSeparation 5 \
#       --chimJunctionOverhangMin 12 \
#       --chimOutType WithinBAM \
#       --chimSegmentReadGapMax 5 \
#       --sjdbGTFfile ${STARREFANNOT} \
#       --twopassMode Basic
#       # --waspOutputMode SAMtag \

module purge
ml picard

java -Xmx16G -jar $EBROOTPICARD/picard.jar SortSam \
      I=${SAMPLENAME}_Aligned.out.bam \
      O=${SAMPLENAME}_qnamesort.bam \
      SORT_ORDER=queryname

java -Xmx16G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=${SAMPLENAME}_qnamesort.bam \
      O=${SAMPLENAME}_qnamesort_markdup.bam \
      M=${SAMPLENAME}_marked_dup_metrics.txt \
      ASSUME_SORT_ORDER=queryname

java -Xmx16G -jar $EBROOTPICARD/picard.jar SortSam \
      I=${SAMPLENAME}_qnamesort_markdup.bam \
      O=${SAMPLENAME}_clean.bam \
      SORT_ORDER=coordinate

module purge
ml SAMtools
samtools index ${SAMPLENAME}_clean.bam


echo "Cleaning up"

module purge
# rm ${SAMPLENAME}*.fastq.gz ${SAMPLENAME}_Aligned.out.bam ${SAMPLENAME}_qnamesort.bam ${SAMPLENAME}_qnamesort_markdup.bam

cd $RNADIR
