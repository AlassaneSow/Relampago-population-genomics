#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=estimate_projections
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=4:00:00
#SBATCH --output=estimate_projections_%j.out
pwd; hostname; date

VCF="/blue/matthewsmith/a.sow/Sebastian_Data/snps/filtered_snps/LD_PRUNED_BIALLELIC_SNPS.vcf.gz"
OUT="/blue/matthewsmith/a.sow/Sebastian_Data/demographic_modelling/"
POP="/blue/matthewsmith/a.sow/Sebastian_Data/snps/filtered_snps/samples_pop.txt"
ml easysfs
cd ${OUT}
easysfs -i ${VCF} -p ${POP} -a -f --preview
