# Variant calling using the Genome Analysis ToolKit 
This code was largely based on the methods described by [Treindl et al. 2023](https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2023.1129867/full), [Tremble et al. 2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521), and notes provided by Dmytro Kryvokhyzha in [GATK: the best practice for genotype calling in a non-model organism](https://evodify.com/gatk-in-non-model-organism/)

## Variant Calling 
### First we call variants using ```HaplotypeCaller```    

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

Then we ran ```HaplotypeCaller``` on all the alignments. See **HaplotypeCaller.sh** 
```console
gatk HaplotypeCaller \   
-R ${REF} \  
-I ${BAM} \    
-ERC GVCF \  
-O $OUT/${SAMPLE}.g.vcf.gz
```    
Then we combinened the gVCF files from HaplotypeCaller into one VCF using ```GenomicsDBImport```. See **gVCF_to_vcf.sh**
```console
gatk GenomicsDBImport \
-R /path_to_reference \
--genomicsdb-workspace-path cohort_db \
--sample-name-map gvcf_map.txt \
```
To make the gvcf_map.txt, simply run this in the HaplotypeCaller output directory  
```console
for f in *.g.vcf.gz; do 
    sample=$(basename "$f" .g.vcf.gz)
    echo -e "${sample}\t$(pwd)/${f}"
done > gvcf_map.txt
```
Index the reference genome and get genomic intervals
```console
gatk CreateSequenceDictionary -R /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna

samtools faidx /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna
awk '{print $1 ":" 1 "-" $2}' /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna.fai > /blue/plp6235/asow/Fomotopsis/ref/genome_intervals.list
```

## Joint-genotyping
### After calling variants, we preform joint genotyping using ```GenotypeGVCFs```. See **joint_genotype.sh**
```console
gatk GenotypeGVCFs \
-R reference.fasta \
-V gendb://cohort_db \
--include-non-variant-sites
-O cohort.vcf.gz
```
## SNP Filtering
### Next, we analyze the quality of the SNPs before filtering. See **snp_extract_and_qc.sh**, **plot_quality.sh**, and **plot_quality.R**
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
-F CHROM -F POS -F QUAL -F QD -GF DP -F MQ -F MQRankSum -F FS -F ReadPOSRankSum -F SOR \
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

### Now, based on the information above and GATK recommendations, we can use ```VariantFiltration```, ```vcftools```, and  ```bcftools``` to remove low quality SNPs. See **filter_seperate_filter_combine.sh**
Please note that because ```pixy``` requires an invariant+variant sites file AND some of the SNP filtering commands (e.g., MQ and QD) remove invariant sites, we would have to filter each site type sepratley. Also the GATK reccommended filters are rather strict and given that the quality of the SNPs are reasonably good, we chose to filter only based on missingness and depth. 

### Missingness and depth filters
```console
vcftools --gzvcf my_vcf.vcf.gz \
--remove-indels \
--max-missing 0.95 \
--min-meanDP 10 \
--max-meanDP 75 \
--recode --stdout | bgzip -c > missingness_filtered_ALLSITES.vcf.gz
```
### Extracting **BIALLELIC SNPs** 
We extract just the SNPs using vcftools
```console
gatk SelectVariants \
  -R ${REF} \
  -V ${OUT}/CLEAN_ALLSITES.vcf.gz \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  -O ${OUT}/RAW_BIALLELIC_SNPS.vcf.gz
```
### LD pruning biallelic SNPs
```console
plink \
--vcf BIALLELIC_SNPS.vcf.gz \
--double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 \
--out ${OUT}/PRUNED_BIALLELIC_SNPS.vcf.gz
```
We started with X biallelic SNPs
```console
bcftools view -H RAW_BIALLELIC_SNPS.vcf.gz | wc -l
```
X sites are included in our PIXY analysis
```console
bcftools view -H pixy_input_CLEAN_ALLSITES.vcf.gz | wc -l
```
and there are X LD pruned bi-allelic SNPs
```console
bcftools view -H PRUNED_BIALLELIC_SNPS.vcf.gz | wc -l
```

## Now we have a file containing bi-allelic SNPs that we can use for ADMIXTURE and PCA and a file with invariant+variant sites for PIXY. With this data we first [linkage-disequillibrium pruned](/Linkage%20Disequilibrium/README.md) the SNPs before analyzing the [population strucutre](/Population%20Structure/Population%20Structure.md) of _Parvodontia relampaga_. 


