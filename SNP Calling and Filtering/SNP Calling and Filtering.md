# Variant calling using the Genome Analysis ToolKit 

## Required Files
Input - BAM file containing reads   
Reference - Reference genome 
## Required Programs*
Genome Analysis Toolkit ```gatk```
```console
module load gatk/4.6.2.0  
```
VCFtools ```vcftools```
```console
module load vcftools/0.1.16
```
BCFtools ```bcftools```
```console
module load bcftools/1.22
```
*If you are curious about the syntax for any of these of programs or want to understand how they work, you can visit the GATK website as they explain it much better than I could (: 
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

## After calling variants, we need to preform joint genotyping using ```GenotypeGVCFs``` 
```console
gatk GenotypeGVCFs \
  -R reference.fasta \
  -V gendb://cohort_db \
  -O cohort.vcf.gz
```
## And extract just the SNPs using ```SelectVariants```
```console
gatk SelectVariants \
-R reference.fasta \
-V cohort.vcf.gz \
--select-type-to-include SNP \
-O unfiltered_SNPS.vcf.gz
```
## Prior to any analysis we need to filter SNPs using ```VariantFiltration```, ```vcftools```, and  ```bcftools```

 We use ```VariantFiltration``` to remove SNPs with low quality metrics according to the GATK Best Practices Workflow

```console
SNP ="path to variant output"
REF ="/path_to_reference"  
```

```console
gatk VariantFiltration \
-R ${REF} \
-V ${SNP} \
-O ${SNP}/filtered_SNPS.vcf \
--filter-name "MQ" - filter "MQ < 40" \
--filter-name "FS" --filter-expression "FS > 60" \
--filter-name "QD" --filter-expression "QD < 2"
```

Afterwards we use ```vcftools``` and ```bcftoools``` to remove uninformative SNPs (E.G., singletons or those with excessive coverage etc)

```console
remove those with >5% missing datat, fewer than two allesl, and read depth below 10 and above 75
```

## Now that we've made a file containing informative, high-quality SNPs we can start analyzing the data. I chose to analyze the [population strucutre](/Population%20Structure/Population%20Structure.md) first. 


