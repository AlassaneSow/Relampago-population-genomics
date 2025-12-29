# Genome assembly 
50X coverage short read whole genome sequences. 
use aaftf 
## Required Files
## Trimmimng and quality control
We trimmed reads with ```trimmomatic v.40```  
``
``
## Mapping to refernce genome
We mapped reads to the P. relampaga reference genome using ```bwa 0.7.19```  
``
``
We then converted the sam files to bams using ```sambamba```  
``
``
We assesed mapping quality and filtered reads using ```sambamba```  
``
``
Mapping statistics after fultering 
``
``
Estimate mean coverage of each individual using ```bedtools```  
