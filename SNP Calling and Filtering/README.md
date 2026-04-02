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
We also have to index the reference genome like so
```console
samtools faidx /path/to/refernce/reference.fna
```
```console
gatk CreateSequenceDictionary \
   -R /path/to/refernce/reference.fna \
   -O /path/to/refernce/reference.dict
```

Then we ran ```HaplotypeCaller``` on all the alignments. Each variable is explained below. 
```console
REF="/path/to/reference.fna"
names="/path/to/names.txt"
OUT="path/to/snp_calling_folder"
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${names}")
BAM="/path/to/filtered_assemblies"
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
## Next, we analyze the quality of the SNPs before filtering. 
To do this we first extract just SNPs from the cohort.vcf.gz
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

## Now, based on the information above and GATK recommendations, we can use ```VariantFiltration```, ```vcftools```, and  ```bcftools``` to remove low quality SNPs.   
Please note that we use two different filtering schemes here. We **1) used site-level and genotype-level filters to make a file that contains just SNPs** (this will be used for most of our analyses) and **2) used less strict genotype-level filters that make sure we don't remove any invariant sites** (this output will be used for pixy which requires a file with variant & invariant sites).

### Filtering for **SNPs ONLY**
We use ```VariantFiltration``` to remove low quality SNPs
```console
SNP ="path to variant output"
REF ="/path_to_reference"  
```
```console
gatk VariantFiltration \
-R ${REF} \
-V ${SNP} \
-O ${SNP}/all_sites_labled.vcf.gz \
--filter "MQ < 40.0" --filter-name "MQ" \
--filter "FS > 60.0" --filter-name "FS" \
--filter "QD < 2" --filter-name "QD" \
```
```console
bcftools view -f PASS all_sites_labeld.vcf.gz -Oz -o cohort_all_sites_VFfiltered.vcf.gz
```
Afterwards we use ```vcftools``` and ```bcftoools``` to remove uninformative SNPs. We...
* kept only chromosonal SNPs
* removed those with >5% missing data
* removed those with read depth <10
* removed those with read depth >75
* Isolates with >80% missing data

To remove non-chromosomal sites we run
```console
bcftools view \
-r Chr1,Chr2..... \
filtered_SNPS.vcf.gz -Oz -o nuclear_filtered_SNPS.vcf.gz
```
Note: In the command above you can use -R and give a list of chromosome names rather than typing them out.

To filtering based on quality we run
```console
vcftools \
--gzvcf nuclear_filtered_SNPS.vcf.gz \ 
--mac 1 \ #this ensures we only keep SNPs
--min-meanDP 10 \
--max-missing 0.70 \
--max-meanDP 75 \
--min-alleles 2 \
--max-alleles 2 \
--recode --recode-INFO-all \
--stdout | bgzip > complete_filtered_SNPS.vcf.gz
```
We started with X SNPs
```console
grep -vc "^#" cohort_all_SNPs.vcf.gz
```
and after all is said and done we are left with X SNPs
```console
grep -vc "^#" complete_filtered_SNPs.vcf.gz
```

### Filtering for **invariant sites AND SNPs**
To calculate the population statistics in ```pixy``` we need to keep both SNPs and invariant sites. Unfortunatley, the filtering above can remove invariant sites. So, for the pixy input, we filter just invariant sites and then combine that output with the file we made above.

```console
vcftools \
--gzvcf .vcf.gz \
--remove-indels \
--max-maf 0 \ #this ensures we keep only invariant sites
--minGQ 20 \
--minDP 5 \
--min-meanDP 8 \
--max-meanDP 80 \
--max-missing 0.80 \
--recode --recode-INFO-all \
--stdout | bgzip > complete_filtered_SNPS.vcf.gz
```
vcftools --gzvcf my_vcf.vcf.gz \
--remove-indels \
--max-missing 0.8 \
--min-meanDP 20 \
--max-meanDP 500 \
--recode --stdout | gzip -c > my_filtered_vcf.vcf.gz

```console
bcftools view \
-r Chr1,Chr2..... \
filtered_SNPS.vcf.gz -Oz -o nuclear_filtered_SNPS.vcf.gz
```
# COMBINE THIS INVARINT FILE WITH THE SNP FILE
## Now that we've made files containing informative, high-quality SNPs we can start analyzing the data. I chose to [LD prune](/Linkage%20Disequilibrium/README.md) the data first before analyzing the [population strucutre](/Population%20Structure/Population%20Structure.md). 



NOTE THIS MIGHT NOT BE NECESSARY! 
First we make a new file that **contains only SNPs** using ```bcftools```
```console
bcftools view -v snps cohort_all_sites_PASS.vcf.gz \
  -Oz -o cohort_variants_only.vcf.gz
```

