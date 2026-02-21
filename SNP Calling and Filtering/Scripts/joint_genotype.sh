#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --job-name=joint_genotype
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=joint_genotype_%j.out
pwd; hostname; date

REF="/path_to_reference"
cohort="path/to/cohort_db"
OUTDIR="path/to/joint/genotype/output"
mkdir -p $OUTDIR
module load gatk

gatk GenotypeGVCFs \
  -R ${REF} \
  -V gendb://${cohort} \
  -O ${OUTDIR}/cohort.vcf.gz
