#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=Hcaller
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=Hcaller_%j.out
pwd; hostname; date

REF="/path_to_reference"
BAM_LIST="/path_to_bam_list.txt"
OUT="/path_to_output_gvcf_folder"
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")
SAMPLE=$(basename "$BAM" .bam)

module load gatk

gatk HaplotypeCaller \
  -R "${REF}" \
  -I "${BAM}" \
  -ERC GVCF \
  -O "${OUT}/${SAMPLE}.g.vcf.gz"
