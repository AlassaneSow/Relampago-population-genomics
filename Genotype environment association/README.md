https://popgen.nescent.org/2018-03-27_RDA_GEA.html

RDA and Univariate associations also annotate loci. Based on Sonstebo 2022.

This is largely based on the Detecting multilocus adaptation using Redundancy Analysis (RDA) vignette in [Population Genetics in R](https://popgen.nescent.org/2018-03-27_RDA_GEA.html). 

First we load the necessary packages
```r
library(psych)
library(vegan)
library(adegenet)
```
We read in the LD-pruned dataset using `adegenet`
```r
data <- read.PLINK("path/to/LD/pruned/data")
```
There can be no NAs in the dataframe so we check and remove them like so
```r
sum(is.na(data))
```
Next we load in the environmental variables and crop them so we just include those in the known range
```r
tiffiles <- list.files("S:/Smith Lab/Alassane/Main/Relampago Blight/WorldClim/wc2.1_2.5m_bio/", pattern = ".tif$", full.names = TRUE)
bioclim <-rast(tiffiles)
se_extent <- ext(-91,-75,25,37)
bioclim_se <- crop(bioclim, se_extent)
```
Then we remove variables that are highly correlated
```r
set.seed(12)
sample_bioclim <- spatSample(bioclim_se, size=10000,
                             method= "random", 
                             na.rm = TRUE, 
                             as.df = TRUE)
cor_matrix <- cor(sample_bioclim, method = "pearson")
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.cex = 0.7, tl.col = "black")
highCorVars <- findCorrelation(cor_matrix, cutoff = 0.7)
bioclim_reduced <- bioclim_se[[-highCorVars]]
removed <- as.data.frame(names(bioclim_se)[highCorVars]) 
names(bioclim_reduced)
writeRaster(bioclim_reduced, 
            filename = "bioclim_reduced.tif", 
            overwrite = TRUE)
```
These variables are...
TABLE
Nowe we can run the RDA using `vegan`
```r

```
