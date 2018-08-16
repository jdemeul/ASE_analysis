#!/bin/bash
#
#SBATCH --job-name=aml_htseq
#SBATCH --output=aml_htseq_%A_%a.out
#SBATCH --error=aml_htseq_%A_%a.err
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6G
#
#SBATCH --cpus-per-task=2
#SBATCH --array=0-14

SAMPLES=( "37" "38" "21" "25" "22" "40" "23" "27" "28" "29" "31" "33" "34" "41" "35" )

/home/jdemeul/projects/mansour/aml_rnaseq_htseq.sh ${SAMPLES[$SLURM_ARRAY_TASK_ID]}
