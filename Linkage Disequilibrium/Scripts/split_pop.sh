#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=split_pop
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu
#SBATCH --ntasks=1
#SBATCH --array=1-#of populations here
#SBATCH --cpus-per-task=4
#SBATCH --mem=2gb
#SBATCH --time=1:00:00
#SBATCH --output=split_pop_%j.out
pwd; hostname; date

module load vcftools#version number

VCF="path/to/vcf/file"
POP="path/to/pop"
P=${POP[$SLURM_ARRAY_TASK_ID]}
NAME=$(basename $P .pop)

vcftools --vcf $VCF \
         --keep $POP \
         --recode --recode-INFO-all \
         --out $NAME
