# Genome assembly 
50X coverage short read whole genome sequences. 
use aaftf 
## Required Files  
All of the variables used here (i.e., {data}) are written out and explained in the full scripts.  
## Quality control and trimming  
Check read quality with ```fastqc```
```console
module load fastqc
fastqc path/to/fastq.gz -o path/to/output/
```
We trimmed reads with ```trimmomatic v.40```. Read quality will determine the trimming options used below. However, assuming everything looks acceptable you can remove adaptors, sequences with <3 base quality at the begining or end of the read, and sequences <50 bp.  
```console
###note to self that {data} refers to the file path and {name} refers to the unique name of each collection
trimmomatic PE -threads 2 -phred33 ${}

trimmomatic PE -threads 2 -phred33 ${data}20180730.A-AT_pool1_${name}_R1.fastq.gz

${data}20180730.A-AT_pool1_${name}_R2.fastq.gz ${out}/${name}.R1.trim.fq.gz ${out}/logs/${name}.R1.un.fq.gz ${out}/${name}.R2.trim.fq.gz ${out}/logs/${name}.R2.un.fq.gz ILLUMINACLIP:${adaptor}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 > ${out}/logs/${name}.trimmo.log  2> ${out}/logs/${name}.trimmo.err

```  
## Mapping to refernce genome
We mapped reads to the P. relampaga reference genome using ```bwa 0.7.19```  
```
```  
We then converted the sam files to bams using ```sambamba```  
```
```  
We assesed mapping quality and filtered reads using ```sambamba```  
```
```  
Mapping statistics after filtering   
```
```  
Estimate mean coverage of each individual using ```bedtools```  
```
```
