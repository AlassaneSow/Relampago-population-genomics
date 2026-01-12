# SNP calling using GATK 

## Required Files
Input - BAM file containing reads 
Reference - Reference genome 

## First need to call variants using GATK ```HaplotypeCaller```    
Load GATK  
``  
module load gatk  
``  
Run ```HaplotypeCaller```  
``
gatk Haplotypecaller \  
-R /path_to_reference \  
-L /path_to_reads  
``    
Then we need to combine the gVCF files from HaplotypeCaller into one VCF using ```CombineGVCFs``` 

Next we need to joint call the SNPs using ```GenotypeGVCFs``` 

Use ```VariantFiltration``` default settings to filter variants

Split output by sample

Select only SNPs from output

