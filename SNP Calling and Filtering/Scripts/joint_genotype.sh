#!/bin/bash
#SBATCH --account=plp6235
#SBATCH --qos=plp6235-b
#SBATCH --job-name=joint_genotype
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=joint_genotype_%j.out
pwd; hostname; date

cd 
ml gatk

gatk GenotypeGVCFs \
-R /blue/plp6235/asow/ref/reference.fna \
-V gendb://cohort_db \
--include-non-variant-sites \
-O /blue/plp6235/asow/snp_calling/snps/raw_genotypes.vcf.gz
