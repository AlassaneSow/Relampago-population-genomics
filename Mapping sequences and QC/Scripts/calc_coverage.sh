#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --job-name=calc_coverage
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=3gb
#SBATCH --time=12:00:00
#SBATCH --output=calc_coverage_%j.out
set -euo pipefail
pwd; hostname; date

ml bedtools
names="path_to_names_file"
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${names}")
quality=20
data="path/to/filtered_assemblies"
OUTDIR="/path/to/coverage"

bedtools genomecov \
-ibam ${data}/${name}_sort_${quality}_dup.bam \
-d | awk '{ total += $3 } END { print total/NR }' \
> ${OUTDIR}/${name}_coverage

cd ${OUTDIR}
echo -e "SampleID\tCoverage" > combined_coverage.tsv
for f in *_coverage; do
  sample=$(echo $f | sed 's/_coverage//')   
  coverage=$(cat "$f")                      
  echo -e "${sample}\t${coverage}" >> combined_coverage.tsv;
done
