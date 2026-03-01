## Phylogenetic trees (Chow)  
### Required files 
FASTA alignment with name up to 256 charcters long. To align WGS with SNP data, use https://github.com/broadinstitute/broad-fungalgroup

Follow these from Tremble et al. Coalescent phylogenomic reconstruction of the group from 702 single-copy gene trees. 

Construct ML phylogenies using ```RAxML v8.0``` using GTRCAT nucleotide substitution model and boostrap based on 1,000 replicates. 

## Maximum clade credibility phylogenetic tree using BEAST (Chow)  
First we need to estimate a mutation rate for each clade using TempEst
``  
``
Then using BEAST under strict molecular clock and GTR substitution model run for 500 million MCMC steps we can make tree
## Evoloutionary rate and molecular dating of clades (Chow)

## Cross-coalescence analysis to time differntiation between clades (Schweizer)

## Time to last common ancestor calculated in BEAST (Chow)

## Coalescent simulations (Stonesenbo)

## Neighbor joining trees with APE and ADEGENET
Neighbor joining trees are important becasue...and they tell you....
To use our SNP in adegenet, we first have to convert the LD pruned data into a genlight object.
Important note from the authors: This function requires
the data to be saved in PLINK using the ‘-recodeA‘ option (see details section
in ?read.PLINK)
```r
library(adegenet)
data <- read.PLINK("path/to/LD/pruned/data")
```
We can use `seploc` to create a list of smaller genlight objects and `dist` to compute parwise distances if there are any issues with memory allocation on the HPC like so:
```r
data_list <- seploc(data, n.block=10, parallel=FALSE)
distances <- lapply(data_list, function(e) bitwise.dist(e))
D <- Reduce("+",distances)
```
Otherwise we simply coompute pairwise distnaces of the entire dataset using `dist`
```r
library(poppr)
D <- bitwise.dist(data)
```
We can now use this to make a NJ tree with `ape`. We colored the tips based on the results of the PCA analysis
```r
library(ape)
library(ggtree)
library(tidyverse)

clusters <- pcafile$cluster
tip_data <- data.frame(label = names(clusters), cluster = factor(clusters))
tree <- nj(D)

ggtree(tree) %<+% tip_data +
geom_tippoint(aes(color=cluster), size = 5)+
theme_tree2()+
scale_color_manual(values=c("red","white"))
```


