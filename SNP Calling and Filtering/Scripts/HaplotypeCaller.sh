#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=Hcaller
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=6:00:00
#SBATCH --output=Hcaller_%j.out
pwd; hostname; date

REF="/path/to/reference.fna"
names="/path/to/names.txt"
OUT="path/to/snp_calling_folder"
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${names}")
BAM="/path/to/filtered_assemblies"

module load gatk

gatk HaplotypeCaller \
  -R ${REF} \
  -I ${BAM}/${name}_sort_20_dup.bam \
  -ERC GVCF \
  -O "${OUT}/${name}.g.vcf.gz"
