#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=map_reads
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=map_reads_%j.out
pwd; hostname; date

names="path_to_names_file"
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${names}")
ref="path/to/reference.fasta"
data="path/to/raw/reads/folder"
OUTDIR="path/to/output/folder"
TMPDIR="${SLURM_TMPDIR}/${name}"
mkdir -p "$TMPDIR"

mkdir -p "${OUTDIR}"
mkdir -p "${TMPDIR}"

module load bwa
module load sambamba

bwa mem -M -t ${SLURM_CPUS_PER_TASK} \
  -R "@RG\tID:${name}\tSM:${name}\tPL:Illumina" \
  ${ref} \
  ${data}/${name}_1.trim.fq.gz \
  ${data}/${name}_2.trim.fq.gz \
| sambamba view -t ${SLURM_CPUS_PER_TASK} -S -f bam -l 0 /dev/stdin \
| sambamba sort /dev/stdin \
    -t ${SLURM_CPUS_PER_TASK} \
    -l 0 \
    -m 10GB \
    --tmpdir ${TMPDIR} \
    -o "${OUTDIR}/${name}.sorted.bam"
