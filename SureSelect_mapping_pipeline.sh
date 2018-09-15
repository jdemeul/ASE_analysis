#!/bin/bash

SAMPLENAME=$1

SUREDIR=/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/SureSelect/Sureselect_FAT1_CAGE_MiSeq_Run-87135524/FASTQ_Generation_2018-08-09_04_06_57Z-114617590

BWAINDEX=/srv/shared/vanloo/pipeline-files/human/references/alignment/GRCh37-1000G/bwa_index/genome.fa
REFGENOME=/srv/shared/vanloo/pipeline-files/human/references/alignment/GRCh37-1000G/genome.fa
TARGETPANEL=/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/SureSelect/Sureselect_FAT1_CAGE_MiSeq_Run-87135524/Sureselect_FAT1_CAGE_MiSeq_Run-87135524_target.bed
KNOWNINDELS=/srv/shared/vanloo/pipeline-files/human/references/b37indels_and_snps/Mills_and_1000G_gold_standard.indels.b37.vcf.bgz
KNOWNSNPS=/srv/shared/vanloo/pipeline-files/human/references/b37indels_and_snps/gnomad_data/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr4.vcf.bgz

echo "Running on sample ${SAMPLENAME}"

FWDREADS=$(find ${SUREDIR} -name "${SAMPLENAME}*_L001_R1_001.fastq.gz")
REVREADS=$(find ${SUREDIR} -name "${SAMPLENAME}*_L001_R2_001.fastq.gz")

CWD=$(pwd)

OUTDIR=/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/SureSelect/Sureselect_FAT1_CAGE_MiSeq_Run-87135524/mapped/${SAMPLENAME}_mapped

mkdir -p ${OUTDIR}
cd ${OUTDIR}


echo "Running trimmomatic"
module purge
ml Trimmomatic

java -Xmx16G -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.38.jar \
      PE -threads 8 -phred33 \
      ${FWDREADS} ${REVREADS} \
      ${SAMPLENAME}_R1_trimmed.fastq.gz ${SAMPLENAME}_R1_trimmed_unpaired.fastq.gz \
      ${SAMPLENAME}_R3_trimmed.fastq.gz ${SAMPLENAME}_R3_trimmed_unpaired.fastq.gz \
      ILLUMINACLIP:${EBROOTTRIMMOMATIC}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


echo "Running bwa-mem"
ml BWA
bwa mem -M -a -t 8 -R "@RG\tID:${SAMPLENAME}\tPL:ILLUMINA\tLB:${SAMPLENAME}\tSM:${SAMPLENAME}" \
      ${BWAINDEX} \
      ${SAMPLENAME}_R1_trimmed.fastq.gz ${SAMPLENAME}_R3_trimmed.fastq.gz \
      > ${SAMPLENAME}.sam


echo "Running samtools/picard"
module purge
ml SAMtools
ml picard

samtools sort -@ 4 -m 4G -o ${SAMPLENAME}_sorted.bam ${SAMPLENAME}.sam

samtools index -@ 4 ${SAMPLENAME}_sorted.bam

java -Xmx16G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=${SAMPLENAME}_sorted.bam \
      O=${SAMPLENAME}_sorted_unique.bam \
      M=${SAMPLENAME}_marked_dup_metrics.txt \
      ASSUME_SORTED=true

samtools index -@ 4 ${SAMPLENAME}_sorted_unique.bam


echo "Running GATK base quality recalibration"
module purge
ml GATK
ml picard

gatk --java-options "-Xmx32G" BaseRecalibrator \
      -I ${SAMPLENAME}_sorted_unique.bam \
      -R ${REFGENOME} \
      --known-sites ${KNOWNSNPS} \
      -L ${TARGETPANEL} \
      -O recal_data.table \
      --interval-padding 100

gatk --java-options "-Xmx32G" ApplyBQSR \
      -R ${REFGENOME} \
      -I ${SAMPLENAME}_sorted_unique.bam \
      --bqsr-recal-file recal_data.table \
      -O ${SAMPLENAME}_sorted_unique_recal.bam


echo "Running Mutect2 and filters"
gatk --java-options "-Xmx32G" Mutect2 \
      -R ${REFGENOME} \
      -I ${SAMPLENAME}_sorted_unique_recal.bam \
      -tumor ${SAMPLENAME} \
      --germline-resource ${KNOWNSNPS} \
      --af-of-alleles-not-in-resource 0.000025 \
      --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
      -L ${TARGETPANEL} \
      --interval-padding 100 \
      -O ${SAMPLENAME}_unmatched_m2_snvs_indels.vcf.gz \
      -bamout ${SAMPLENAME}_clean_m2realign.bam
      #   --panel-of-normals ${BASEDIR}/8samplePoN.vcf.gz \

# still need FilterMutectCalls
gatk --java-options "-Xmx16G" FilterMutectCalls \
      -V ${SAMPLENAME}_unmatched_m2_snvs_indels.vcf.gz \
      --max-germline-posterior .9 \
      -O ${SAMPLENAME}_unmatched_m2_snvs_indels_filt.vcf.gz

java -Xmx4G -jar $EBROOTPICARD/picard.jar BedToIntervalList \
      I=${TARGETPANEL} \
      O="${TARGETPANEL/.bed/.list}" \
      SD="${REFGENOME/.fa/.dict}"

# this command requires columns 4 and 5 of BED file to be filled out (or bed transformed into list using picard's BedToIntervalList)
gatk --java-options "-Xmx8G" CollectSequencingArtifactMetrics \
      -I ${SAMPLENAME}_sorted_unique_recal.bam \
      -O ${SAMPLENAME}_SequencingArtifactMetrics \
      -R ${REFGENOME} \
      --DB_SNP ${KNOWNSNPS} \
      --INTERVALS ${TARGETPANEL/.bed/.list}

# and FilterByOrientationBias / Oxo-G etc
gatk --java-options "-Xmx8G" FilterByOrientationBias \
      --artifact-modes 'G/T' \
      -V ${SAMPLENAME}_unmatched_m2_snvs_indels_filt.vcf.gz \
      -P ${SAMPLENAME}_SequencingArtifactMetrics.pre_adapter_detail_metrics \
      -O ${SAMPLENAME}_unmatched_m2_snvs_indels_filt_x2.vcf.gz

gatk --java-options "-Xmx4G" SelectVariants \
    -R ${REFGENOME} \
    -V ${SAMPLENAME}_unmatched_m2_snvs_indels_filt_x2.vcf.gz \
    -O ${SAMPLENAME}_unmatched_m2_snvs_indels_filt_x2_pass.vcf.gz \
    -L ${TARGETPANEL} \
    --exclude-filtered true

# make SV/indel calls using SvABA
module purge
ml SvABA

svaba run -p 4 --germline \
      -t ${SAMPLENAME}_sorted_unique_recal.bam \
      -k ${TARGETPANEL} \
      -a ${SAMPLENAME} \
      -G ${BWAINDEX} \
      --dbsnp-vcf ${KNOWNSNPS}

echo "Cleaning up"

module purge
rm ${SAMPLENAME}*_trimmed.fastq.gz ${SAMPLENAME}*_trimmed_unpaired.fastq.gz \
    ${SAMPLENAME}.sam ${SAMPLENAME}_sorted.bam* \
    ${SAMPLENAME}_sorted_unique.bam* ${SAMPLENAME}_marked_dup_metrics.txt

cd ${CWD}
