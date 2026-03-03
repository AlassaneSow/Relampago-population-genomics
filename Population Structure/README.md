# PCA and STRUCTURE 
normal stuff and   
LD prune the dataset using PLINK.

study effect of geographic distance on populations
## Required Files
## Required Programs

To determine the population strucutre of _P. relampaga_ we used principal component analyses (PCA) with ```ADGENT``` and Bayesian clustering (STRUCTURE) with ```LEA```
## Principal component analysis (PCA)
First we load the LD-pruned data as a genlight object. PLINK needs to have recodeA
```r
library(adegenet)
data <- read.PLINK("path/to/LD/pruned/data")
```
Then we run the PCA
```r
PCA <- glPca(data)

ggplot(PCA, aes(X=PCA1, Y=PCA2, fill=population))+
geom_point(shape=21, size=3)+
theme_bw
```
Eigenvalues graph below.
We can compare the observed realtionship to the NJ tree. See here for more infomration. 
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

## STRUCTURE
First we use the `snmf2` function 
```r
snmf_result <- snmf(data, K=1:10, ploidy=2, entropy=T, alpha=100,project="new")
```

admixture proportions on map
https://bookdown.org/hhwagner1/LandGenCourse_book/WE_9.html
