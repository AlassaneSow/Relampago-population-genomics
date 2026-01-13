# SNP calling using GATK 

## Required Files
Input - BAM file containing reads   
Reference - Reference genome 

## First we need to call variants using GATK ```HaplotypeCaller```    
Load GATK  
```  
module load gatk  
```
Because ```HaplotypeCaller``` can only handle one sample at a time we need to loop the command for each alignment. To do this we create a list file containing every alignment.   
```
ls /path_to_alignments/*bam > bam_list.txt
```

Then we run ```HaplotypeCaller``` on all the alignments
```
REF=/path_to_reference  
BAM_LIST=/path_to_bam_list.txt  
OUT=/path_to_output_gvcf folder
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BAM_LIST)
SAMPLE=$(basename "$BAM" .bam)
gatk HaplotypeCaller \   
-R $REF \  
-I $BAM \    
-ERC GVCF \  
-O $OUT/${SAMPLE}.g.vcf.gz
```    
Then we need to combine the gVCF files from HaplotypeCaller into one VCF using ```GenomicsDBImport```

```
gatk GenomicsDBImport \
-R ref

```
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

