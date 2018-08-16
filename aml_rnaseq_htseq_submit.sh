#!/bin/bash
#
#SBATCH --job-name=aml_htseq
#SBATCH --output=aml_htseq_%A_%a.out
#SBATCH --error=aml_htseq_%A_%a.err
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#
#SBATCH --cpus-per-task=2
#SBATCH --array=0-6

# SAMPLES=( "37" "38" "21" "25" "22" "40" "23" "27" "28" "29" "31" "33" "34" "41" "35" )
SAMPLES=( "HL60" "MOLM16" "OCI_AML2" "OCI_AML3" "TF1" "THP1_S6" "THP1_S13" )

/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/aml_rnaseq_htseq.sh ${SAMPLES[$SLURM_ARRAY_TASK_ID]}
