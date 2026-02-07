#!/bin/bash
#SBATCH --array=1-100

REF="/path_to_reference"
BAM_LIST="/path_to_bam_list.txt"
OUT="/path_to_output_gvcf_folder"

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")
SAMPLE=$(basename "$BAM" .bam)

gatk HaplotypeCaller \
  -R "${REF}" \
  -I "${BAM}" \
  -ERC GVCF \
  -O "${OUT}/${SAMPLE}.g.vcf.gz"
