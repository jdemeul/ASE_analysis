#!/bin/bash
#
#SBATCH --job-name=aml_map_cell
#SBATCH --output=aml_map_cell_%A_%a.out
#SBATCH --error=aml_map_cell_%A_%a.err
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6G
#
#SBATCH --cpus-per-task=10
#SBATCH --array=0-5

# SAMPLES=( "HL60" "MOLM16" "OCI_AML2" "OCI_AML3" "TF1" "THP1_S6" "THP1_S13" )
SAMPLES=( "MOLM16" "OCI_AML2" "OCI_AML3" "TF1" "THP1_S6" "THP1_S13" )

/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/eml_aml_cell_lines_RNASeq.sh ${SAMPLES[$SLURM_ARRAY_TASK_ID]}
