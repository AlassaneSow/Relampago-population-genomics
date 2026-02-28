# Calculating population summary statistics and identifyinng local adaptation
This code was largely based on the methods described by [Tremble et al. 2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521). Here we calculate nucleotide diversity (&pi;), the fixation index (F<sub>ST</sub>), absolute genetic divergence (d<sub>XY</sub>), and Taijama's D. To calculate sliding window size, see https://speciationgenomics.github.io/ld_decay/
## Required files
Populations file
* The populactions file is necessary for fst and dxy.  
variant sites VCF  
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
The populations file is necessary to calculate F<sub>ST</sub> and d<sub>XY</sub>. To make this file we used the populations defined by the PCA. Below is an example of the file.
```
ColS1  SF
Fak1  SF
NATL4  NF
```
## Identifying highly divergent loci and genes
To identify highly divergent loci, we calculated F<sub>ST</sub> and d<sub>XY</sub> for each gene in each lineage. We also used PCadapt which identifies locally adapted loci.   
### F<sub>ST</sub> and d<sub>XY</sub>
First we calculated F<sub>ST</sub> and d<sub>XY</sub> only in genic regions of the genome. To do this we created a .bed file containing the cordinates of all of the genic regions by filtering the GTF file to only include genes  
```console
awk '$3=="gene"{print $1"\t"$4"\t"$5}' path/to/reference.gff > genes.bed
```  
We calculated the statistics as follows  
```console
pixy --stats fst dxy \  
--vcf path/tp/data.vcf.gz \  
--populations PCA_populations.txt \  
--bed_file genes.bed \  
```
To ensure we only retain significant divergent loci we only kept loci in the top 5% of both F<sub>ST</sub> and d<sub>XY</sub>
```R
# To run this see pixy_filter.R and use the Rscript command
library(tidyverse)

fst <- read_tsv("path/to/pixy_fst_output.tsv")
dxy <- read_tsv("path/to/pixy_dxy_output.tsv")
genes <- read_tsv("path/to/genes.bed",
                  col_names = c("chromosome","start","end","gene_id"))
data <- fst %>%
  select(chromosome,start,end,pop1,pop2,avg_wc_fst) %>%
  left_join(
  dxy %>%
    select(chromsome, start,end,pop1,pop2,avg_dxy),
  by = c("chromsome", "start","end","pop1","pop2","avg_dxy")
)
data <- data %>%
  left_join(genes,
            by = c("chromosome","start","end"))
divergent_loci <- data %>%
  group_by(pop1,pop2) %>%
  filter(
    avg_wc_fst >= quantile(avg_wc_fst,0.95,na.rm=TRUE),
    avg_dxy >= quantile(avg_dxy, 0.95, na.rm=TRUE)) %>%
    ungroup()
write_tsv(data, "raw_divergent_loci.tsv")
write_tsv(divergent_loci, "divergent_loci_fst_dxy.tsv")
```
### pcadapt
First we have to create a .vcf that contains just SNPs in genic regions. 
```console
bedtools intersect -a input.vcf -b genes.gff -header -wa > output.vcf
```
Then we have to convert the .vcf to .bed so pcadapt can use it. 
```console
plink \
--vcf ${VCF} \
--make-bed \
--out ${OUT}/genotype_pcadapt
```
Load in the data
```R
library(pcadapt)
library(ggplot2)
library(tidyverse)
path <- "path_to_directory/genotype_pcadapt.bed"
data <- read.pcadapt(path, type = "bed")
```
Prior to running pcadapt, we need to define the number of ```K``` principal components. We will use Cattell's rule to find the number of PCs that correspond to population strucure. For more information on this [see here](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html)
```R
K_graph <- pcadapt(input=data, K=20)
plot(K_graph, option="screeplot")
```
This shows that we should retain X PCs
PICTURE HERE
Now we can run pcadapt with the appropiate number of PCs
```R
locally_adapted <- pcadapt(data, K=X)
```
We checked for the frequency of SNPs with significant p-values
```R
ggplot(locally_adapted, aes(pvalues))+
  geom_histogram(bins = 50, color = "purple")+
  geom_vline(xintercept = 0.05, color = "black")
```
We used Bonferroni correction to retain only SNPs that are reliably and significanly locally adapted. We used alpha=0.01
```R
sig_locally_adapted <- locally_adapted %>%
mutate(p_adjusted = p.adjust(pvalues, method="bonferroni")) %>%
filter(p_adjusted > 0.01)
write_tsv(sig_locally_adapted, "pcadapt_output.tsv")
```

