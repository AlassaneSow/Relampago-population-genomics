# Genome assembly 
We used COMPANY Illlumina HiSeq4000 platform to sequence 50X coverage whole genomes. This code is largley based on the work in [Treindl et al. 2023](https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2023.1129867/full)
### Required Files  
Reference genome from JGI 
Raw sequence data
### Required programs
FASTQC ```fastqc```   
Trimmomatic ```trimmomatic v.40```  
Burrows-Wheeler Aligner ```bwa 0.7.19```  
Sambamba ```sambamba```   
Bedtools ```bedtools```  
## First we check the sequencing quality using ```fastqc```  
```console
fastqc path/to/fastq.gz -o path/to/output/
```
It looks like the quality is good. See NAME.html for more read quality information. 
PUT IMAGE HERE!

## Next we trimmed reads with ```trimmomatic v.40```.  
The read quality will  will determine which trimming options one uses. Because our quality looks okay we simply removed...
* adaptors
* sequences with <3 base quality at the begining or end of the read
* sequences <50 bp  

```console
data="path/to/raw/reads/folder/"
name="path_to_names_file.txt"
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
```
### Note  
To make the names file, we ran this in the folder containing either the forward or reverse reads
```console
ls sharedprefix_*_R1.fastq.gz | \
  sed 's/^sharedprefix_//' | \
  sed 's/_R1.fastq.gz$//' > names.txt
```
## Next we mapped reads to the _P. relampaga_ reference genome using ```bwa 0.7.19```  
```
bwa index path/to/reference.fa
```
```
bwa mem -M -t 2 \
${Ref} \
${data}${name}_1.trim.fq.gz \
${data}${name}_2.trim.fq.gz \
-M \
-R "@RG\tID:${name}\tSM:${name}\tPL:Illumina" > ${OUTDIR}/${name}.sam
```  
We then converted the sam files to bams using ```sambamba```  
```
sambamba view -t 4 \
-S ${TMPDIR}/${name}.sam \
-f bam \
-l 0 \
-o /dev/stdout|sambamba sort\
/dev/stdin \
-o /dev/stdout \
-t ${proc} \
-l 0 \
-m 10GB \
--tmpdir ${TMPDIR} > ${OUTDIR}/${name}.sort.bam
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
