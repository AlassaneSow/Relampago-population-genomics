https://popgen.nescent.org/2018-03-27_RDA_GEA.html
https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html

RDA and Univariate associations also annotate loci. Based on Sonstebo 2022.

This is largely based on the Detecting multilocus adaptation using Redundancy Analysis (RDA) vignette in [Population Genetics in R](https://popgen.nescent.org/2018-03-27_RDA_GEA.html). 
Before using R, we have to reformat the data so it can be used. In the terminal we run: 
```console
plink --tped yournameoffile.tped --tfam yournameoffile.tfam --recode A --out yournameoffile
```
First we load the necessary packages
```r
library(psych)
library(vegan)
library(adegenet)
```
We read in the LD-pruned dataset and remove NA values 
```r
data <- read.table("path/to/LD/pruned/data.raw", header=TRUE)
data <- data[,6:37029] #remove metadata coloumns

sum(is.na(data))# X NAs in the matrix (~X% missing data)

data_NA <- data #copy df
data_NA[is.na(data_NA)] <- colMeans(data_NA, na.rm=TRUE)[col(data_NA)][is.na(data_NA)] #replace missing SNPs with mean genotype of its coloumn 
genotypes <- data_NA
sum(is.na(genotypes)) #no more missing data
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
We standardize the envrionmental data 
```r
bioclim_reduced <- scale(bioclim_reduced)
```
For better readability later I suggest changing the variable names here like so
```r
DO ME
```
Nowe we can run the RDA using `vegan` DOUBLE CHECK HOW TO DO THIS WITH RASTER DATA
```r
RDA <- rda(data ~., data=as.data.frame.(bioclim_reduced))
```
We calculate R-squared and other summary statistics as follows
R-squared tells you the % of genomic variation explained by the constained ordinations variable 
```r
RsquareAdj(RDA) #R^2=0.00X, meaning the ordination explains 0.0X% of the variation in the data 
```
The eigenvalues for the axes tells us about the variance explained by each axis
```r
summary(eigenvals(RDA, model="constrained"))
```
We can also check the significance using an anova test
```r
significance_full <- anova.cca(RDA, parallel=getOption("mc.cores'))
significance_per_axis <- anova.cca(RDA, by=axis, parallel=getOption("mc.cores'))
```
We can quickly make a plot like so
```r
plot(RDA, scaling=3)
```
The red points are all the SNPs, the grey points are each individual, and the vectors are the environmental predictors.  
Now we can plot just the SNPs like so
```r
plot(RDA,type="n", scaling=3)
points(RDA, display="species", pch=20, cex=0.7, col="gray32", scaling=3)#the "species" option plots just SNPs
text(RDA, scaling=3, display="bp", col="#0868ac", cex=1)                          
```
and just the individuals like so
```r
plot(RDA,type="n", scaling=3)
points(RDA, display="sites", pch=20, cex=0.7, col="gray32", scaling=3)#the "species" option plots just SNPs
text(RDA, scaling=3, display="bp", col="#0868ac", cex=1) 
```
To make these plots a bit more informative we can use the code below to color each individual by its designated population.
```r
###Color individuals by population
populations <- read.table("populations.txt", header=TRUE) #load populations file we made earlier
site_scores <- scores(RDA, display="sites") #extract the site scores
df <- as.data.frame(site_scores)
df$sample <- rownames(df)
RDA_ggplot <- merge(df, populations, by="sample", all.x=TRUE) #merge RDA results with population names
```
and replot like this
```r
ggplot(data=rda)
```
Individual SNPs linked to these envriomental variables are..
```
loadings <- scores(rda_model, choices=1:2, display="species")

Calculate SD threshold:

z <- apply(loadings,2,scale)

outliers <- which(abs(z) > 3, arr.ind=TRUE)

Or per axis:

sd1 <- sd(loadings[,1])
mean1 <- mean(loadings[,1])

outliers_axis1 <- which(abs(loadings[,1] - mean1) > 3*sd1)
```
