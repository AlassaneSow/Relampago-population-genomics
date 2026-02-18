#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=vcf_to_plink
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu
#SBATCH --ntasks=1
#SBATCH --array=1-#number of populations -1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4gb
#SBATCH --time=2:00:00
#SBATCH --output=vcf_to_plink_%j.out
pwd; hostname; date

module load plink

VCFS=(path/to/folder/*.recode.vcf)
VCF=${VCFS[$SLURM_ARRAY_TASK_ID]}
POP=$(basename $VCF .recode.vcf)

plink --vcf $VCF \
      --allow-extra-chr \
      --keep-allele-order \
      --make-bed \
      --out ${POP}
