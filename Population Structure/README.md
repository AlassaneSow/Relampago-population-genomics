# PCA, STRUCTURE, isolation by distance 
## Required Files
## Required packages
```r
library(adegenet)
library(vcfR)
library(tidyverse)
library(caret)
library(vegan)
library(LEA)
```
## PCA
To determine the population strucutre of _P. relampaga_ we used principal component analyses (PCA) with ```ADGENT``` and Bayesian clustering (STRUCTURE) with ```LEA```
## Principal component analysis (PCA)
# NOTE TO SELF TO ADD SEBASTIAN ET AL WORKBOOK NOTES TO THE SCRIPTS FOLDER
First we load the LD-pruned data and convert to a genlight object. PLINK needs to have recodeA
```r
raw_data <- read.vcfR(path/to/vcf)
data <- vcfR2genlight(raw_data)
```
We then remove the file paths from the sample names
```r
indNames(data) <- gsub("/media/hdd2tb/popgen/Pilze/PopGen/BAM-Ac/|\\.bam|'|/media/hdd2tb/popgen/Pilze/PopGen/BAM-Fp/|\\.bam|","",indNames(data))
```
Next we run the PCA. Note that if you run this in a slurm script you MUST add the nf option as `adegenet` will prompt you to choose how many principal components to report. If you do not do this the job will never finish.
```r
PCA <- glPca(data,parallel = TRUE,nf=10)

ggplot(PCA, aes(X=PCA1, Y=PCA2, fill=population))+
geom_point(shape=21, size=3)
```
INCLUDE PCA IMAGE AND EIGNENVALUES PLOTs.  

We can compare the observed realtionships to those in the NJ tree. To see how this figure was made, see here. 
```r
library(ape)
library(poppr)
tree_data <- bitwise.dist(data)
tre <- nj(tree_data)

ggtree(tre) %<+% tip_data +
geom_tippoint(aes(color=cluster), size = 5)+
theme_tree2()+
scale_color_manual(values=c("red","white"))
```
To see how we made the final PCA for the publication, see PCA.R
## STRUCTURE
First we load the vcf, extract geontypes, and use the vcf_to_lea function to convert to a format LEA can use. 
```r
vcf <- read.vcfR(".vcf.gz")
genotypes <- extract.gt(vcf, element="GT", as.numeric=TRUE)

vcf_to_lea <- function(genotypes, output_file) {
  # Convert missing data (NA) to 9 (LEA's missing data code)
  genotypes[is.na(genotypes)] <- 9
  # Write to a .geno file
  write.table(genotypes, file = output_file, sep = "", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
vcf_to_lea(geontypes, "LEA_genotypes_for_strucutre.geno")
obs <- snmf(".geno", K = 1:10, repetitions = 10, entropy = TRUE, project = "new", CPU = 12)
```
Note that you can always re-load the snmf output and that you should only have to run the analysis once. 
```
obs <- load.snmf(".smnfproject)
```
admixture proportions on map
https://bookdown.org/hhwagner1/LandGenCourse_book/WE_9.html

## IBD
calcuclate identity by state disrtacne matirxi between all pairs of undividuales uinsg plink with the --distance 1-ibs command. See Krah et al. 
Also calcualte pairwise genome-wide Fst and plot against distance as in Treindl et al.
