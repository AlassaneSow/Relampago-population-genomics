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

Then we run ```HaplotypeCaller``` on all the alignments. See HCall.sh for the full script. 
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
-R /path_to_reference \
--genomicsdb-workspace-path cohort_db \
--sample-name-map gvcf_map.txt \
```
To make the gvcf_map.txt, simply run this in the HaplotypeCaller output directory  
```
for gvcf in *.g.vcf.gz; do
  sample=$(basename "$gvcf" .g.vcf.gz)
  echo -e "${sample}\t$(pwd)/${gvcf}"
done > gvcf_map.txt
```  

Next we need to joint call the SNPs using ```GenotypeGVCFs``` 
```
gatk GenotypeGVCFs \
  -R reference.fasta \
  -V gendb://cohort_db \
  -O cohort.vcf.gz
```

Use ```VariantFiltration``` default settings to filter variants

Split output by sample

Select only SNPs from output

