# Variant calling using the Genome Analysis ToolKit 
This code was largely based on the methods described by [Treindl et al. 2023](https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2023.1129867/full), [Tremble et al. 2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521), and notes provided by Dmytro Kryvokhyzha in [GATK: the best practice for genotype calling in a non-model organism](https://evodify.com/gatk-in-non-model-organism/)

### Required Files
Input - BAM file containing reads   
Reference - Reference genome 
### Required Scripts
HCall.sh
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
--include-non-variant-sites
-O cohort.vcf.gz
```
## Prior to any analysis we need to filter SNPs using ```VariantFiltration```, ```vcftools```, and  ```bcftools```. 

### Prior to any filtering we checked the quality of the SNPs. First we extract just SNPs from the vcf.gz
```console
gatk SelectVariants \
-R reference.fasta \
-V cohort.vcf.gz \
-selectType SNP \
-o cohort_all_SNPs.vcf.gz
```
Then we extract the quality metrics from the vcf.gz
```console
gatk VariantsToTable \
-R reference.fasta \
-V cohort_all_SNPs.vcf.gz \
-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPOSRankSum -F SOR \
--allowMissingData \
-o SNP_Scores.table
```
Lastly, we plot the quality metrics in R.
```r
library('gridExtra')
library('ggplot2')

SNPs <- read.csv("path/to/SNP_Scores.table")

DP <- ggplot(SNPS, aes(x=DP, fill=purple))+
      geom_density(alpha=0.3)+
      geom_vline(xintercept=c(10,6200))

QD <- ggplot(SNPS, aes(x=QD, fill=purple))+
      geom_density(alpha=0.3)+
      geom_vline(xintercept=2, size=0.7)

FS <- ggplot(SNPS, aes(x=FS, fill=purple))+
      geom_density(alpha=0.3)
      geom_vline(xintercept=60, size=0.7)

MQ <- ggplot(SNPS, aes(x=MQ, fill=purple))+
      geom_density(alpha=0.3)+
      geom_vline(xintercept=40, size=0.7)

MQRankSum <- ggplot(SNPS, aes(x=MQRankSum, fill=purple))+
      geom_density(alpha=0.3)

SOR <- ggplot(SNPS, aes(x=SOR, fill=purple))+
      geom_density(alpha=0.3)
      geom_vline(xintercept=4, size=0.7)

ReadPosRankSum <- ggplot(SNPS, aes(x=ReadPosRankSum, fill=purple))+
      geom_density(alpha=0.3)+
      geom_vline(xintercept=c(10,-10), size=0.7)

svg("path/to/output/QualityPlots.svg", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, nrow=4)
dev.off()
```

### We used ```VariantFiltration``` to remove SNPs with low quality metrics as defined in the GATK Best Practices Workflow. 

```console
SNP ="path to variant output"
REF ="/path_to_reference"  
```
```console
gatk VariantFiltration \
-R ${REF} \
-V ${SNP} \
-O ${SNP}/all_sites_filtered.vcf.gz \
--filter "MQ < 40.0" --filter-name "MQ" \
--filter "FS > 60.0" --filter-name "FS" \
--filter "QD < 2" --filter-name "QD" \
```
```console
bcftools view -f PASS all_sites_filtered.vcf.gz -Oz -o cohort_all_sites_PASS.vcf.gz
```
Afterwards we use ```vcftools``` and ```bcftoools``` to remove uninformative SNPs. We...
* kept only chromosonal SNPs
* removed those with >5% missing data
* removed those with read depth <10
* removed those with ead depth >75
* Isolates with >80% missing data

Remove non-chromosomal sites
```console
bcftools view \
-r Chr1,Chr2..... \
filtered_SNPS.vcf.gz -Oz -o nuclear_filtered_SNPS.vcf.gz
```
Note: In the command above you can use -R and give a list of chromosome names rather than typing them out.

Filtering based on quality
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
For most of the analyses (e.g., PCA and ADMIXTURE) we need  a file that contains only SNPs. We extracted them from the vcf.gz using ```bcftools```
```console
bcftools view -v snps cohort_all_sites_PASS.vcf.gz \
  -Oz -o cohort_variants_only.vcf.gz
```
To calculate the population statistics in ```pixy``` we need to keep both SNPs and invariant sites so we don't do anything. 

We started with X SNPs
```console
grep -vc "^#" cohort_all_SNPs.vcf.gz
```
and after all is said and done we are left with X SNPs
```console
grep -vc "^#" complete_filtered_SNPs.vcf.gz
```

## Now that we've made file containing informative, high-quality SNPs we can start analyzing the data. I chose to [LD prune](/Linkage%20Disequilibrium/README.md) the data first before analyzing the [population strucutre](/Population%20Structure/Population%20Structure.md). 


