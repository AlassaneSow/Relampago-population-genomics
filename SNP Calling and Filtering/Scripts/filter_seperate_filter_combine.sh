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
#SBATCH --output=extract_and_qc_%j.out
pwd; hostname; date

REF="/blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna"
RAW="/blue/matthewsmith/a.sow/Sebastian_Data/snps/raw_genotypes.vcf.gz"
OUT="/blue/matthewsmith/a.sow/Sebastian_Data/snps/filtered_snps"

ml gatk
ml samtools
ml vcftools
ml bcftools 

## Filter based on depth and missingness #max missing can min DP can be relaxed if needed.
vcftools --gzvcf ${RAW} \
--max-missing 0.95 \
--min-meanDP 10 \
--max-meanDP 75 \
--recode --stdout | bgzip -c > ${OUT}/CLEAN_ALLSITES.vcf.gz

tabix -p vcf ${OUT}/CLEAN_ALLSITES.vcf.gz

## Extract ONLY BIALLELIC SNPs
gatk SelectVariants \
  -R ${REF} \
  -V ${RAW} \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  -O ${OUT}/RAW_BIALLELIC_SNPS.vcf.gz

##Index files, combine, and sort
cd ${OUT}
#tabix -p vcf CLEAN_ALLSITES.vcf.gz
#tabix -p vcf RAW_BIALLELIC_SNPS.vcf.gz

bcftools sort CLEAN_ALLSITES.vcf.gz -Oz -o pixy_input_CLEAN_ALLSITES.vcf.gz
bcftools index pixy_input_CLEAN_ALLSITES.vcf.gz
X=$(bcftools view -H pixy_input_CLEAN_ALLSITES.vcf.gz | wc -l)
echo "There are $X sites in this file"

bcftools sort RAW_BIALLELIC_SNPS.vcf.gz -Oz -o BIALLELIC_SNPS.vcf.gz
bcftools index BIALLELIC_SNPS.vcf.gz
Y=$(bcftools view -H BIALLELIC_SNPS.vcf.gz | wc -l)
echo "There are $Y sites in this file"

