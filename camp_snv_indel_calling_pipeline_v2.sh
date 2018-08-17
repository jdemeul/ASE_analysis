#!/bin/bash

SAMPLENAME=$1
#RUNDIR=$2
BASEDIR=/camp/lab/vanloop/working/demeulj/projects/2016_mansour_ASE_T-ALL/data/SureSelect2017_all_data_files/SureSelect-20Sep17-48078054/
#FWDREADS=$3
#REVREADS=$4

RUNDIR=$(find ${BASEDIR} -type d -name "${SAMPLENAME}*")
FWDREADS=$(find ${RUNDIR} -name "${SAMPLENAME}*_L001_R1_001.fastq.gz")
REVREADS=$(find ${RUNDIR} -name "${SAMPLENAME}*_L001_R2_001.fastq.gz")

BWAREFGENOME=/camp/lab/vanloop/working/demeulj/reference/hg38/Homo_sapiens_assembly38.fasta
PICARD=/home/camp/demeulj/.local/easybuild/software/picard/2.17.8-Java-1.8.0_131/picard.jar
GATK=/home/camp/demeulj/gatk-4.0.1.1/gatk


echo "Running on sample ${SAMPLENAME}"
cd ${RUNDIR}


# echo "Running trimmomatic"
# module purge
# ml Trimmomatic

# java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 4 -phred33 \
# ${FWDREADS} ${REVREADS} \
# ${SAMPLENAME}_R1_trimmed.fastq.gz ${SAMPLENAME}_R1_trimmed_unpaired.fastq.gz \
# ${SAMPLENAME}_R2_trimmed.fastq.gz ${SAMPLENAME}_R2_trimmed_unpaired.fastq.gz \
# ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# echo "Running bwa-mem"
# ml BWA
# bwa mem -t 4 ${BWAREFGENOME} ${SAMPLENAME}_R1_trimmed.fastq.gz ${SAMPLENAME}_R2_trimmed.fastq.gz > ${SAMPLENAME}.sam


# echo "Running samtools/picard"
# module purge
# ml SAMtools
ml Java

# samtools view -@ 4 -b ${SAMPLENAME}.sam -o ${SAMPLENAME}.bam

# java -jar ${PICARD} SortSam \
#       I=${SAMPLENAME}.bam \
#       O=${SAMPLENAME}_qnamesort.bam \
#       SORT_ORDER=queryname

# java -jar ${PICARD} MarkDuplicates \
#       I=${SAMPLENAME}_qnamesort.bam \
#       O=${SAMPLENAME}_qnamesort_markdup.bam \
#       M=${SAMPLENAME}_marked_dup_metrics.txt \
#       ASSUME_SORT_ORDER=queryname

# java -jar ${PICARD} AddOrReplaceReadGroups \
#       I=${SAMPLENAME}_qnamesort_markdup.bam \
#       O=${SAMPLENAME}_clean.bam \
#       SORT_ORDER=coordinate \
#       RGID=1 \
#       RGLB=lib1 \
#       RGPL=illumina \
#       RGPU=unit1 \
#       RGSM=${SAMPLENAME}

# samtools index ${SAMPLENAME}_clean.bam

# echo "Running Mutect2 and filters"
# ${GATK} --java-options "-Xmx8G" Mutect2 \
#         -R ${BWAREFGENOME} \
#         -I ${SAMPLENAME}_clean.bam \
#         -tumor ${SAMPLENAME} \
#         --germline-resource /home/camp/demeulj/af-only-gnomad.hg38.vcf.gz \
#         --panel-of-normals ${BASEDIR}/8samplePoN.vcf.gz \
#         --af-of-alleles-not-in-resource 0.0000025 \
#         --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
#         -L /camp/lab/vanloop/reference/Genomics/gatk/hg38/wgs_calling_regions.hg38.interval_list \
#         -O ${SAMPLENAME}_unmatched_m2_snvs_indels.vcf.gz \
#         -bamout ${SAMPLENAME}_clean_m2realign.bam

## only for PoN
## /home/camp/demeulj/gatk-4.0.1.1/gatk CreateSomaticPanelOfNormals \
##     -vcfs AA14185440-57724764/AA14185440_unmatched_m2_snvs_indels.vcf.gz \
##     -vcfs DW11338750-57725778/DW11338750_unmatched_m2_snvs_indels.vcf.gz \
##     -vcfs GK11338175-57725779/GK11338175_unmatched_m2_snvs_indels.vcf.gz \
##     -vcfs ME1333016-57737706/ME1333016_unmatched_m2_snvs_indels.vcf.gz \
##     -vcfs JET11333536-57736713/JET11333536_unmatched_m2_snvs_indels.vcf.gz \
##     -vcfs KB14185227-57721789/KB14185227_unmatched_m2_snvs_indels.vcf.gz \
##     -vcfs MH11337537-57721791/MH11337537_unmatched_m2_snvs_indels.vcf.gz \
##     -vcfs TP14185528-57720801/TP14185528_unmatched_m2_snvs_indels.vcf.gz \
##     -O 8samplePoN.vcf.gz


# # still need FilterMutectCalls
# ${GATK} --java-options "-Xmx8G" FilterMutectCalls \
#         -V ${SAMPLENAME}_unmatched_m2_snvs_indels.vcf.gz \
#         --max-germline-posterior .9 \
#         -O ${SAMPLENAME}_unmatched_m2_snvs_indels_filt.vcf.gz

# # and FilterByOrientationBias
# ${GATK} --java-options "-Xmx8G" CollectSequencingArtifactMetrics \
#       -I ${SAMPLENAME}_clean.bam \
#       -O ${SAMPLENAME}_artifact \
#       --FILE_EXTENSION ".txt" \
#       -R ${BWAREFGENOME}

# ${GATK} --java-options "-Xmx8G" FilterByOrientationBias \
#       -AM G/T \
#       -V ${SAMPLENAME}_unmatched_m2_snvs_indels_filt.vcf.gz \
#       -P ${SAMPLENAME}_artifact.pre_adapter_detail_metrics.txt \
#       -O ${SAMPLENAME}_unmatched_m2_snvs_indels_filt_x2.vcf.gz

${GATK} --java-options "-Xmx4G" SelectVariants \
    -R ${BWAREFGENOME} \
    -V ${SAMPLENAME}_unmatched_m2_snvs_indels_filt_x2.vcf.gz \
    -O ${SAMPLENAME}_unmatched_m2_snvs_indels_pass.vcf.gz \
    -L chr5:149909849-150359849 \
    --exclude-filtered true

echo "Cleaning up"

module purge
# rm ${SAMPLENAME}*_trimmed.fastq.gz ${SAMPLENAME}*_trimmed_unpaired.fastq.gz \
#     ${SAMPLENAME}.sam ${SAMPLENAME}.bam ${SAMPLENAME}_qnamesort.bam \
#     ${SAMPLENAME}_qnamesort_markdup.bam ${SAMPLENAME}_unmatched_m2_snvs_indels_filt.vcf.gz* 
# tar -czf ${SAMPLENAME}_artifact_metrics.tar.gz *_metrics.txt
# zcat ${SAMPLENAME}_unmatched_m2_snvs_indels_filt_x2.vcf.gz | grep -P "^#|PASS" | gzip > ${SAMPLENAME}_unmatched_m2_snvs_indels_filt_x2_clean.vcf.gz
cd $BASEDIR
