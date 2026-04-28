#!/bin/bash
#SBATCH --account=plp6235
#SBATCH --qos=plp6235-b
#SBATCH --job-name=extract_and_qc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50gb
#SBATCH --time=8:00:00
#SBATCH --output=extract_and_qc_%j.out
pwd; hostname; date

ml gatk 

gatk SelectVariants \
-R /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna \
-V /blue/plp6235/asow/Fomotopsis/snp_calling/snps/raw_genotypes.vcf.gz \
-select-type SNP \
-O /blue/plp6235/asow/Fomotopsis/snp_calling/snps/raw_snps.vcf.gz

gatk VariantsToTable \
-R /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna \
-V /blue/plp6235/asow/Fomotopsis/snp_calling/snps/raw_snps.vcf.gz \
-F CHROM -F POS -F QUAL -F QD -GF DP -F MQ -F MQRankSum -F FS -F ReadPOSRankSum -F SOR \
-O /blue/plp6235/asow/Fomotopsis/snp_calling/snps/raw_snp_scores.table
