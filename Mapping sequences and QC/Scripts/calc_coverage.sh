#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --job-name=calc_coverage
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=calc_coverage_%j.out
set -euo pipefail
pwd; hostname; date
module load bedtools
names="path_to_names_file"
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${names}")
quality=20
OUTDIR="path/to/folder/with/${name}_sort_${quality}_dup.bam/"
mkdir -p ${OUTDIR}stats!${quality}

bedtools genomecov \
-ibm ${OUTDIR}${name}_sort_${quality}_dup.bam \
-d | awk'{ total += $3 } END { print total/NR }' \
> ${OUTDIR}stats!${quality}/${name}_cov
