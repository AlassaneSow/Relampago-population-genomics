# PCA and STRUCTURE 
normal stuff and   
LD prune the dataset using PLINK.

study effect of geographic distance on populations
## Required Files
## Required Programs

To determine the population strucutre of _P. relampaga_ we used principal component analyses (PCA) with ```ADGENT``` and Bayesian clustering (STRUCTURE) with ```ADMIXTURE```
## Principal component analysis (PCA)
First we load the LD-pruned data as a genlight object. PLINK needs to have recodeA
```r
library(adegenet)
data <- read.PLINK("path/to/LD/pruned/data")
```
Then we run the PCA
```r
PCA <- glPca(data)
```

## STRUCTURE
