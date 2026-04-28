#!/bin/bash
#SBATCH --account=plp6235
#SBATCH --qos=plp6235-b
#SBATCH --job-name=filter_seperate_allsites_snps
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50gb
#SBATCH --time=8:00:00
#SBATCH --output=extract_and_qc_%j.out
pwd; hostname; date

ml gatk
ml samtools
ml vcftools
ml bcftools 

##Broad filters to remove sites with extremley high & low depth and missingness. #need to see if 75  reasonable for invariant sites in my data
vcftools --gzvcf my_vcf.vcf.gz \
--remove-indels \
--max-missing 0.95 \
--min-meanDP 10 \
--max-meanDP 75 \
--recode --stdout | bgzip -c > missingness_filtered_ALLSITES.vcf.gz

##Extract ONLY SNPs and filter them based on quality
vcftools --gzvcf missingness_filtered_ALLSITES.vcf.gz \
--mac 1 \
--recode --stdout | bgzip -c > missingness_filtered_SNPS.vcf.gz

gatk VariantFiltration \
-R ${REF} \
-V missingness_filtered_SNPS.vcf.gz \
-O quality_labeled_SNPS.vcf.gz \
--filter "MQ < 40.0" --filter-name "MQ" \
--filter "FS > 60.0" --filter-name "FS" \
--filter "QD < 2" --filter-name "QD"

bcftools view -f PASS quality_labeled_SNPS.vcf.gz -Oz -o quality_filtered_SNPS.vcf.gz

##Extract ONLY BI-ALLELIC SNPs
vcftools --gzvcf quality_filtered_SNPS.vcf.gz \
--min-alleles 2 \
--max-alleles 2 \
--recode --stdout | bgzip -c > quality_filtered_BIALLELIC_SNPS.vcf.gz

##Extract ONLY invariant sites
gatk SelectVariants \
  -R ${REF} \
  -V missingness_filtered_ALLSITES.vcf.gz \
  --select-type-to-include NO_VARIATION \
  -O quality_filtered_INVARIANT.vcf.gz

##Index files, combine, and sort
tabix -p vcf quality_filtered_INVARIANT.vcf.gz
tabix -p vcf quality_filtered_SNPS.vcf.gz

bcftools concat \
--allow-overlaps \
quality_filtered_INVARIANT.vcf.gz quality_filtered_SNPS.vcf.gz \
-O z -o quality_filtered_ALLSITES.vcf.gz

bcftools sort quality_filtered_ALLSITES.vcf.gz -Oz -o pixy_input_quality_filtered_ALLSITES.vcf.gz
bcftools index pixy_input_quality_filtered_ALLSITES.vcf.gz
