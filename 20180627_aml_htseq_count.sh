#!/bin/bash

SAMPLES=( "37" "38" "21" "25" "22" "40" "23" "27" "28" "29" "31" "33" "34" "41" "35" "CD34" )

for i in $(seq 1 8);
do
htseq-count --format=bam --order=pos --stranded=reverse --quiet --type=exon --idattr=gene_id /srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/AML/RNA/mapped/${SAMPLES[i]}_mapped/${SAMPLES[i]}_clean.bam /srv/shared/vanloo/pipeline-files/human/references/annotation/GENCODE/gencode.v23lift37.annotation.gtf > /srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_RNAseq_counts/${SAMPLES[i]}.alignments.bam.count &
done
wait

echo "done waiting 1"

for i in $(seq 9 16);
do
htseq-count --format=bam --order=pos --stranded=reverse --quiet --type=exon --idattr=gene_id /srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/AML/RNA/mapped/${SAMPLES[i]}_mapped/${SAMPLES[i]}_clean.bam /srv/shared/vanloo/pipeline-files/human/references/annotation/GENCODE/gencode.v23lift37.annotation.gtf > /srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_RNAseq_counts/${SAMPLES[i]}.alignments.bam.count &
done
wait

echo "done waiting 2"
