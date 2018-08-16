#!/bin/bash

SAMPLENAME=$1

htseq-count --format=bam --order=pos --stranded=reverse --quiet --type=exon --idattr=gene_id /srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/AML/RNA/mapped/${SAMPLENAME}_mapped/${SAMPLENAME}_clean.bam /srv/shared/vanloo/pipeline-files/human/references/annotation/GENCODE/gencode.v23lift37.annotation.gtf > /srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_RNAseq_counts/${SAMPLENAME}.alignments.bam.count
