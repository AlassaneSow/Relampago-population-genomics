# Variant calling using the Genome Analysis ToolKit 

## Required Files
Input - BAM file containing reads   
Reference - Reference genome 
## Required Programs  
Genome Analysis Toolkit ```gatk```
```console
module load gatk/4.6.2.0  
```
## First we need to call variants using ```HaplotypeCaller```    

Because ```HaplotypeCaller``` can only handle one sample at a time we need to loop the command for each alignment. To do this we create a list file containing every alignment.   
```console
ls /path_to_alignments/*bam > bam_list.txt
```

Then we run ```HaplotypeCaller``` on all the alignments. See HCall.sh for the full script. 
```console
REF="/path_to_reference"  
BAM_LIST="/path_to_bam_list.txt"  
OUT="/path_to_output_gvcf folder"
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BAM_LIST)
SAMPLE=$(basename "$BAM" .bam)
```
```console
gatk HaplotypeCaller \   
-R ${REF} \  
-I ${BAM} \    
-ERC GVCF \  
-O $OUT/${SAMPLE}.g.vcf.gz
```    
Then we need to combine the gVCF files from HaplotypeCaller into one VCF using ```GenomicsDBImport```  
```console
gatk GenomicsDBImport \
-R /path_to_reference \
--genomicsdb-workspace-path cohort_db \
--sample-name-map gvcf_map.txt \
```
To make the gvcf_map.txt, simply run this in the HaplotypeCaller output directory  
```console
for gvcf in *.g.vcf.gz; do
  sample=$(basename "$gvcf" .g.vcf.gz)
  echo -e "${sample}\t$(pwd)/${gvcf}"
done > gvcf_map.txt
```  

## After calling variants, we need to joint preform joint genotyping using ```GenotypeGVCFs``` 
```console
gatk GenotypeGVCFs \
  -R reference.fasta \
  -V gendb://cohort_db \
  -O cohort.vcf.gz
```
## Now we can extract just SNPs using ```SelectVariants```
```console
gatk SelectVariants \
-R reference.fasta \
-V cohort.vcf.gz \
--select-type-to-include SNP \
-O unfiltered_SNPS.vcf.gz
```
## Prior to any analysis we need to filter SNPs using ```VariantFiltration``` and ```vcftools```

 We use ```VariantFiltration``` to remove SNPs with low quality metrics according to the GATK Best Practices Workflow

```console
SNP ="path to variant output"
```

```console
gatk VariantFiltration \
-R ${REF} \
-V ${SNP} \
-O ${SNP}/filtered_SNPS.vcf \
--filter-name "XXX" --filter-expression "XXXX" \
--filter-name "XXX" --filter-expression "XXXX" \
--filter-name "XXX" --filter-expression "XXXX" \
--filter-name "XXX" --filter-expression "XXXX" \
--filter-name "XXX" --filter-expression "XXXX" 
```

Afterwards we use ```vcftools``` and ```bcftoools``` to remove uninformative SNPs (i.e., singletons, those with excessive coverage etc)  


