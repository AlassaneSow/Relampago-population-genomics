#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=filter_variants_step1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=filter_variants_step1_%j.out
pwd; hostname; date

SNP ="path to variant output"
REF ="/path_to_reference"  
OUT ="path/to/output/inSNP/folder"

module load gatk
module load bcftools

gatk VariantFiltration \
-R ${REF} \
-V ${SNP} \
-O ${OUT}/all_sites_labeled.vcf.gz \
--filter "MQ < 40.0" --filter-name "MQ" \
--filter "FS > 60.0" --filter-name "FS" \
--filter "QD < 2" --filter-name "QD" 

bcftools view -f PASS ${OUT}/all_sites_labeled.vcf.gz -Oz -o ${OUT}/cohort_all_sites_VFfilterd.vcf.gz

bcftools view \
-r Chr1,Chr2..... \
${OUT}/cohort_all_sites_VFfilterd.vcf.gz -Oz -o ${OUT}/cohort_all_sites_VF_nuclear_filtered.vcf.gz
