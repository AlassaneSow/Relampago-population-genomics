# SNP calling using GATK 

## Required Files
Input - BAM file containing reads   
Reference - Reference genome 

## First need to call variants using GATK ```HaplotypeCaller```    
Load GATK  
```  
module load gatk  
```
Then create a list file containing all of the alignments.   
```
ls /path_to_alignments/*bam > bam_list.txt
```  

Run ```HaplotypeCaller```  
```
gatk Haplotypecaller \   
-R /path_to_reference \  
-L /path_to_mapped_reads \
--sample-ploidy 1 \  
-ERC GVCF \  
-O /path_to_output_gvcf folder
```    
Then we need to combine the gVCF files from HaplotypeCaller into one VCF using ```CombineGVCFs```   
```
gatk CombineGVCFs \
-R /path_to_reference \
$(for f in /path_to_output_gvcf; do echo "--variant $f"; done) \ \ #this loops through all the vcf files and prints their names rather than typing them one by one
-O /path_to_output folder
```

Next we need to joint call the SNPs using ```GenotypeGVCFs``` 

Use ```VariantFiltration``` default settings to filter variants

Split output by sample

Select only SNPs from output

