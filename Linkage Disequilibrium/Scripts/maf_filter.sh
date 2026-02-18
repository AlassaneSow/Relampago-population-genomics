#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=maf_filter
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --array=1-#number of populations -1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4gb
#SBATCH --time=2:00:00
#SBATCH --output=maf_filter_%j.out
pwd; hostname; date

module load plink

BEDS=(path/to/folder/*.bed)
BED=${BEDS[$SLURM_ARRAY_TASK_ID]}
POP=$(basename $BED .bed)

plink \
--bfile path/to/bed/output/${POP} \
--allow-extra-chr \
--maf 0.05 \
--make-bed \
--out ${POP}_filt
