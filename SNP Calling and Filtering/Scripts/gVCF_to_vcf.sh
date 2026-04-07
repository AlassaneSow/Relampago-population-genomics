#!/bin/bash
#SBATCH --account=plp6235
#SBATCH --qos=plp6235-b
#SBATCH --job-name=combine_gvcf
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150gb
#SBATCH --time=48:00:00
#SBATCH --output=combine_gvcf_%j.out
pwd; hostname; date

cd /blue/plp6235/asow/Fomotopsis/snp_calling
ml gatk

gatk GenomicsDBImport \
-R /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna \
--genomicsdb-workspace-path fomotopsis_cohort_db \
-L /blue/plp6235/asow/Fomotopsis/ref/genome_intervals.list \
--merge-input-intervals true \
--sample-name-map gvcf_map.txt
