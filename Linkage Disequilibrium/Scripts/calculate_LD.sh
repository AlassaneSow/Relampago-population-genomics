#!/bin/bash
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=calculate_LD
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --array=1-#number of populations -1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4gb
#SBATCH --time=2:00:00
#SBATCH --output=calculate_LD_%j.out
pwd; hostname; date

module load plink

BEDS=(path/to/filtered/folder/*.bed)
BED=${BEDS[$SLURM_ARRAY_TASK_ID]}
POP=$(basename $BED .bed)

plink \
--file path/to/filtered/output/${POP}_filt \
--allow-extra-chr \
--r2 \
--ld-window-kb 10000 \
--ld-window-r2 0 \
--out ${POP}_filt
