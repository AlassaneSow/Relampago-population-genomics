# Calculating population summary statistics. 
This code was largely based on the methods described by [Tremble et al. 2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521). Here we calculate nucleotide diversity (&pi;), the fixation index (F<sub>ST</sub>), absolute genetic divergence (d<sub>XY</sub>), and Taijama's D 
## Required files
Populations file
* The populactions file is necessary for fst and dxy.  
LD pruned VCF  
non-LD pruned VCF
Reference transcriptome
## Required programs
```pixy```
## Calculating statistics
We used ```PIXY``` to calculate &pi;, F<sub>ST</sub>, d<sub>XY</sub>, and Taijama's D as in 10kb sliding windows. 
```console
pixy --stats pi fst dxy tajima_d \  
--vcf path/to/appropriate/file.vcf.gz \  
--populations PCA_populations.txt \  
--window_size 10000 \  
--n_cores 2  
```
The populations file is necessary to calculate F<sub>ST</sub> and d<sub>XY</sub>. To make this file we defined populations as either the geographic populations OR the populations detected by the PCA. We then compared these values to see if they are similar. 
* If they are different I think I should determine Fst between the major geographic areas.
```
ColS1  SF
Fak1  SF
NATL4  NF
```
## Identifying highly divergent loci and genes
To identify genes causing local adaptation we calculated F<sub>ST</sub> and d<sub>XY</sub> only in genic regions of the genome. To do this we created a .bed file containing the cordinates of all of the genic regions by filtering the GTF file to only include genes  
```console
awk '$3=="gene"{print $1"\t"$4"\t"$5}' path/to/reference.gff > genes.bed
```  
We then calculated the statistics as follows  
```console
pixy --stats fst dxy \  
--vcf path/tp/data.vcf.gz \  
--populations PCA_populations.txt \  
--bed_file genes.bed \  
```
To ensure we only retain significant divergent loci we only kept loci in the top 5% of both F<sub>ST</sub> and d<sub>XY</sub>
```R
# To run this see pixy_filter.R and use the Rscript command
data <- read.table("path/to/pixy_output.tsv", header=TRUE, sep="\t")
fst_cut <- quantile(data$fst, 0.95, na.rm=TRUE)
dxy_cut <- quantile(data$dxy, 0.95, na.rm=TRUE)
out <- d[data$fst >= fst_cut & data$dxy >= dxy_cut, ]
write.table(out, "path/to/pixy/output/5perc_genic.tsv, sep="\t", quote=FALSE, row.names=FALSE)
```

Lastly we remove loci that were divergent due to genetic drift using 
```console
```

