#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=gVCF_to_vcf
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=gVCF_to_vcf_%j.out
pwd; hostname; date

REF="/path_to_reference"
gvcf_map="path/to/gvcf_map.txt"

module load gatk

gatk GenomicsDBImport \
-R ${REF} \
--genomicsdb-workspace-path cohort_db \
--sample-name-map ${gvcf_map} \
