# Variant calling using the Genome Analysis ToolKit 
This code was largely based on the methods described by [Treindl et al. 2023](https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2023.1129867/full), [Tremble et al. 2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521), and notes provided by Dmytro Kryvokhyzha in [GATK: the best practice for genotype calling in a non-model organism](https://evodify.com/gatk-in-non-model-organism/)

### Required Files
Input - BAM file containing reads   
Reference - Reference genome 
### Required Programs
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
## First we call variants using ```HaplotypeCaller```    

Because ```HaplotypeCaller``` can only handle one sample at a time we need to call each alignment one at a time. To do this we created a list file containing every alignment.   
```console
ls /path_to_alignments/*bam > bam_list.txt
```

Then we ran ```HaplotypeCaller``` on all the alignments. See HCall.sh for the full script. Each variable is explained below. 
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
Then we combinened the gVCF files from HaplotypeCaller into one VCF using ```GenomicsDBImport```  
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

## After calling variants, we preform joint genotyping using ```GenotypeGVCFs``` 
```console
gatk GenotypeGVCFs \
  -R reference.fasta \
  -V gendb://cohort_db \
  -O cohort.vcf.gz
```
## For most of the analyses (e.g., PCA and ADMIXTURE) we need to only use SNPs. We extracted just the SNPs using ```SelectVariants```
```console
gatk SelectVariants \
-R reference.fasta \
-V cohort.vcf.gz \
--select-type-to-include SNP \
--exclude-non-variants \
--exclude-filtered \
-O unfiltered_SNPS.vcf.gz
```
## To calculate the population statistics in ```pixy``` we need to keep both SNPs and invariant sites. We simply ran the same command without the exclude-non-variants and --select-type-to-include commands. 
```console
gatk SelectVariants \
-R reference.fasta \
-V cohort.vcf.gz \
--exclude-filtered \
-O unfiltered_SNPS.vcf.gz
```
# Note to self to filter quality first then seperate variant and invariant sites.
## Prior to any analysis we need to filter SNPs using ```VariantFiltration```, ```vcftools```, and  ```bcftools```

We used ```VariantFiltration``` to remove SNPs with low quality metrics as defined in the GATK Best Practices Workflow. 

```console
SNP ="path to variant output"
REF ="/path_to_reference"  
```

```console
gatk VariantFiltration \
-R ${REF} \
-V ${SNP} \
-O ${SNP}/filtered_SNPS.vcf \
--filter "MQ < 40.0" --filter-name "MQ" \
--filter "FS > 60.0" --filter-name "FS" \
--filter "QD < 2" --filter-name "QD" \
```

Afterwards we use ```vcftools``` and ```bcftoools``` to remove uninformative SNPs. We...
* kept only chromosonal SNPs
* removed those with >5% missing data
* removed those with read depth <10
* removed those with ead depth >75
* Isolates with >80% missing data

```console
bcftools view \
-r Chr1,Chr2..... \
filtered_SNPS.vcf.gz -Oz -o nuclear_filtered_SNPS.vcf.gz
```
```console
vcftools \
--gzvcf nuclear_filtered_SNPS.vcf.gz \ 
--mac 3 \
--min-meanDP 10 \
--max-missing 0.70 \
--max-meanDP 75 \
--min-alleles 2 \
--max-alleles 2 \
--recode --recode-INFO-all \
--stdout | bgzip > complete_filtered_SNPS.vcf.gz
```

Optional filtering based on the analysis and quality of the data. See [Treindl et al. 2023](https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2023.1129867/full):
* Remove SNPs with >30% missing in one population (I will have to check if missigness is wildly uneven prior to doing this)
* Remove SNPs with >80% missing data
* --mac 3. This removes rare alleles and is likely optional. I Can retain for some analyses but probably important for selective sweep scans? Idk. 

## Now that we've made a file containing informative, high-quality SNPs we can start analyzing the data. I chose to [LD prune](/Linkage%20Disequilibrium/README.md) the data first before analyzing the [population strucutre](/Population%20Structure/Population%20Structure.md). 


