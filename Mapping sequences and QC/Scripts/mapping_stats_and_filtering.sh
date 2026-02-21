#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --job-name=stats_and_filtering
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=stats_and_filtering_%j.out
set -euo pipefail
pwd; hostname; date
module load sambamba
module load gatk

names="path_to_names_file"
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${names}")
data=path/to/input.bams
quality=20
OUTDIR="path/to/output/folder"
TMPDIR="${SLURM_TMPDIR}/${name}"
mkdir -p "${TMPDIR}"
mkdir -p ${OUTDIR}/stats
mkdir -p ${OUTDIR}/stats_dup
mkdir -p ${OUTDIR}/statsQ${quality}
TMP_BAM=${TMPDIR}/${name}_sort_${quality}.bam
DEDUP_BAM=${OUTDIR}/${name}_sort_${quality}_dup.bam

sambamba flagstat ${data}/${name}.sort.bam \
> ${OUTDIR}/stats/${name}

sambamba view \
-F "mapping_quality >= ${quality}" \
-t ${SLURM_CPUS_PER_TASK} \
-f bam \
-l 0 \
${data}/${name}.sort.bam \
-o ${TMP_BAM}

gatk MarkDuplicates \
--TMP_DIR=${TMPDIR} \
--INPUT=${TMP_BAM} \
--OUTPUT=${DEDUP_BAM} \
--METRICS_FILE=${OUTDIR}/stats_dup/${name} \
--VALIDATION_STRINGENCY=LENIENT \
--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024

sambamba index ${DEDUP_BAM}

sambamba flagstat ${DEDUP_BAM} \
> ${OUTDIR}/statsQ${quality}/${name}
