# Measuring linkage disequilibrium (LD) and linkage pruning SNPs
This code is largley based on the methods described by [Tremble et al. 2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521) and the notes provided by Mark Ravinet & Joana Meier in [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/).   
  
Some analyses, such as PCA and admixture models, asssume that SNPs are independent. As such we need to remove SNPs in LD to make accurate interpretations. Futhermore, while there are many ways to choose how to set the sliding window paramater when calculating population genetic statistics (e.g.,&pi;,F<sub>ST</sub>, d<sub>XY</sub>), most agree that the sliding window should be set to the size greater than LD decay. Below we investigate LD decay across all isolates and within the populations defined in our population structure analyes. 

## Loading data
plink --vcf quality_filtered_BIALLELIC_SNPS.vcf.gz \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--geno 0.1 \
--mind 0.3 \
--make-bed \
--out biallelicSNPs

## Linkage pruning
We used ```plink``` to filter out SNPs in high LD (r2>0.2) by scanning 50kb windows in 10 bp steps.
```console
plink \
--vcf ${VCF} \
--double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 \
--out ${OUT}
```
We used this otput in our PCA analysis. 

## Measuring LD decay
We used ```plink``` to measure LD decay across all samples and within each population. 
To calculate across all samples, we simply ran
```console
plink --vcf quality_filtered_BIALLELIC_SNPS.vcf.gz \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 \
--geno 0.1 \
--mind 0.5 \
-r2 gz \
--ld-window 100 \
--ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed \
--out biallelic_LD
```
Load into R 
```
library(readr)
ld <- read_table("your_file.ld.gz")
```
To calculate LD within populations, we 
First, we made seperate vcf files for each population. To do this we made .pop population files that contained the names of each sample in a population.
```console
FID IID
sample1 sample1
sample2 sample2
sample3 sample3
```
Then we used ```vcftools``` to split the .vcf file. 
```console
vcftools \
--vcf $VCF \
--keep $POP \
--recode --recorde-INFO-all \
--out $NAME 
```
Next we converted the .vcf to a format ```plink``` can use. 
```console
plink \
--vcf $VCF \
--allow-extra-chr \
--make-bed \
--out $POP
```
To reduce noise we removed rare SNPs
```console
plink \
--bfile ${POP} \
--allow-extra-chr \
--maf 0.05 \
--make-bed \
--out ${POP}_filt
```
Lastly we calculated LD in 10kb windows
```console
plink \
--file ${POP}_filt \
--allow-extra-chr \
--r2 \
--ld-window-kb 10000 \
--ld-window-r2 0 \
 --out ${POP}_filt
```
# Figure here!

Optional significance test for wether individuals belong to same population
