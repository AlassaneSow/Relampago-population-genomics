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
## First we checked the sequencing quality using ```fastqc```  
```console
fastqc path/to/fastq.gz -o path/to/output/
```
It looks like the quality is good. See NAME.html for more read quality information. 
PUT IMAGE HERE!

## Next we trimmed reads with ```trimmomatic v.40```.  
The read quality will determine which trimming options one uses. Because our quality looks okay we simply removed...
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
```console
bwa index path/to/reference.fa
```
```console
bwa mem -M -t 2 \
${Ref} \
${data}${name}_1.trim.fq.gz \
${data}${name}_2.trim.fq.gz \
-M \
-R "@RG\tID:${name}\tSM:${name}\tPL:Illumina" > ${OUTDIR}/${name}.sam
```  
We then converted the sam files to bams using ```sambamba```  
```console
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
We then assesed the mapping statisitcs, filtered low quality alignments (those with mapping quality less than 20), and removed duplicates using ```sambamba``` and ```gatk```  
```console
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
```  
Lastly we estimated mean coverage of each individual using ```bedtools```  
```console
bedtools genomecov \
-ibm ${OUTDIR}${name}_sort_${quality}_dup.bam \
-d | awk'{ total += $3 } END { print total/NR }' \
> ${OUTDIR}stats!${quality}/${name}_cov
```
