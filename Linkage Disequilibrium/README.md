USe plink and follow this https://speciationgenomics.github.io/ld_decay/
# Measuring linkage disequilibrium (LD) and linkage pruning SNPs
This code is largley based on the methods described by [Tremble et al. 2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521) and the notes provided by Mark Ravinet & Joana Meier in [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/).   
  
Some analyses, such as PCA and admixture models, asssume that SNPs are independent. As such we need to remove SNPs in LD to make accurate interpretations. Futhermore, while there are many ways to choose how to set the sliding window paramater when calculating population genetic statistics (e.g.,&pi;,F<sub>ST</sub>, d<sub>XY</sub>), most agree that the sliding window should be set to the size where LD decays. Below we investigate LD decay across all isolates and within the populations defined in our population structure analyes.   

All of this is to say that some of the code below has to be done before analyzing population structure and some relies on the population structure results. 

## Required Files
filtered vcf.gz

## Required Programs
```plink version number```
```console
module load plinkvXXX
```
## Linkage pruning

## Measuring LD decay

Optional significance teste for wether individuals belong to same population
