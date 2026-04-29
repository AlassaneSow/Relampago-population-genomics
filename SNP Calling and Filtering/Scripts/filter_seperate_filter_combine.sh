#!/bin/bash
#SBATCH --account=plp6235
#SBATCH --qos=plp6235-b
#SBATCH --job-name=filter_seperate_allsites_snps
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50gb
#SBATCH --time=30:00:00
#SBATCH --output=filter_seperate_allsites_snps_%j.out
pwd; hostname; date

REF="/blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna"
RAW="/blue/matthewsmith/a.sow/Sebastian_Data/snps/raw_genotypes.vcf.gz"
OUT="/blue/matthewsmith/a.sow/Sebastian_Data/snps/filtered_snps"

ml gatk
ml samtools
ml vcftools
ml bcftools 

##Broad filters to remove sites with extremley high & low depth and missingness. I need to see if 75 reasonable for invariant sites in my data.
vcftools --gzvcf ${RAW} \
--remove-indels \
--max-missing 0.95 \
--min-meanDP 10 \
--max-meanDP 75 \
--recode --stdout | bgzip -c > ${OUT}/missingness_filtered_ALLSITES.vcf.gz

##Extract ONLY SNPs and filter them based on quality
vcftools --gzvcf ${OUT}/missingness_filtered_ALLSITES.vcf.gz \
--mac 1 \
--recode --stdout | bgzip -c > ${OUT}/missingness_filtered_SNPS.vcf.gz

gatk VariantFiltration \
-R ${REF} \
-V ${OUT}/missingness_filtered_SNPS.vcf.gz \
-O ${OUT}/quality_labeled_SNPS.vcf.gz \
--filter "MQ < 40.0" --filter-name "MQ" \
--filter "FS > 60.0" --filter-name "FS" \
--filter "QD < 2" --filter-name "QD"

bcftools view -f PASS ${OUT}/quality_labeled_SNPS.vcf.gz -Oz -o ${OUT}/quality_filtered_SNPS.vcf.gz

##Extract ONLY BI-ALLELIC SNPs
vcftools --gzvcf ${OUT}/quality_filtered_SNPS.vcf.gz \
--min-alleles 2 \
--max-alleles 2 \
--recode --stdout | bgzip -c > ${OUT}/quality_filtered_BIALLELIC_SNPS.vcf.gz

##Extract ONLY invariant sites
gatk SelectVariants \
  -R ${REF} \
  -V ${OUT}/missingness_filtered_ALLSITES.vcf.gz \
  --select-type-to-include NO_VARIATION \
  -O ${OUT}/quality_filtered_INVARIANT.vcf.gz

##Index files, combine, and sort
cd ${OUT}
tabix -p vcf quality_filtered_INVARIANT.vcf.gz
tabix -p vcf quality_filtered_SNPS.vcf.gz

bcftools concat \
--allow-overlaps \
quality_filtered_INVARIANT.vcf.gz quality_filtered_SNPS.vcf.gz \
-O z -o quality_filtered_ALLSITES.vcf.gz

bcftools sort quality_filtered_ALLSITES.vcf.gz -Oz -o pixy_input_quality_filtered_ALLSITES.vcf.gz
bcftools index pixy_input_quality_filtered_ALLSITES.vcf.gz
