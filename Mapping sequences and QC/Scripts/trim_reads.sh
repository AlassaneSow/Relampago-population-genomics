#!/bin/sh
#SBATCH --account=plantpath
#SBATCH --qos=plantpath
#SBATCH --job-name=trim_reads
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12gb
#SBATCH --time=12:00:00
#SBATCH --output=trim_reads_%j.out
pwd; hostname; date

data="path/to/raw/reads/folder/"
name="path_to_names_file"
out="path/to/output/folder"
adaptor="path/to/adaptor/file.txt"

trimmomatic PE -phred33 -threads 2 \ 
${data}sharedprefix_${name}_R1.fastq.gz \
${data}sharedprefix_${name}_R2.fastq.gz \
${out}/${name}.R1.trim.fq.gz \
${out}/logs/${name}.R1.un.fq.gz \
${out}/${name}.R2.trim.fq.gz \
${out}/logs/${name}.R2.un.fq.gz \
ILLUMINACLIP:${adaptor}:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:50 > ${out}/logs/${name}.trimmo.log \
2> ${out}/logs/${name}.trimmo.err
