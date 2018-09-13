#!/bin/bash
#
#SBATCH --job-name=aml_cell_wxs
#SBATCH --output=aml_cell_wxs%A_%a.out
#SBATCH --error=aml_cell_wxs%A_%a.err
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#
#SBATCH --array=0-3

# SAMPLES=(AA14185440 DF11338384 DW11338750 GK11338175 KAJ11333588 LB11332785 ME1333016 PD11337536 BW11337822 DW11337794 EP14185811 JET11333536 KB14185227 Loucy MH11337537 TP14185528)
# SAMPLES=(CCRF-CEM_S5 DU-528_S3 JURKAT_S4 KOPT-K1_S1 LOUCY_S7 MOLT-3_S9 P12-ICHIKAWA_S6 PF-382_S2 RPMI-8402_S8)
SAMPLES=( "MOLM-16" "OCI-AML2" "OCI-AML3" "TF-1" )

ml Trimmomatic BWA SAMtools picard GATK/3.8-1-Java-1.8.0_162

/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/WXS_mapping_pipeline.sh ${SAMPLES[$SLURM_ARRAY_TASK_ID]}

module purge