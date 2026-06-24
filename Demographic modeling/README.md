# Demographic modeling using fastsimcoal2
This code is based on the [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/fastsimcoal2/), the [fastsimcoal2 manual](https://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal28.pdf), and methods used by Tremble et al. [2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521).
## Preparing data
Before running ```fastsimcoal2``` we have to estimate a site-frequency spectrum using easy SFS. To do this we will need the LD-pruned vcf and a populations file. We will use the popualtions identified in the K=X [ADMIXTURE plots](https://github.com/AlassaneSow/Relampago-population-genomics/tree/main/Population%20Structure#structure).  
  
Although the vcf contains very high quality SNPs, ```easysfs``` can not use missing data and we may have to downsample our data to retain the highest number of sites. If there is lots of missing data or low covergae, the alternate program ANGS can create an SDS and account for genotype uncertanties.
```console
# estimate the best projections
easSFS -i ${VCF} -p ${POP} -a -f --preview
```
The output looks like this. Look out for where the number of segregating sites sharply decreases, this will be the number we choose for the --proj flag. THe proj flag lists the number of individuals in both populations. Because we filtered for missing data and removed poorly assembles isolates, we do not have to downsample.
Next we created the SFS files
```console
easysfs -i ${VCF} -p ${POP} -a -f --proj 34,34
```
## Modeling
We chose to model X scenarios to explain the observed P. relampaga population structure (see below). The .tpl file we used for each scenatio is in Files/Inputs/ and the paramters are explained below. 

### Two populations with no migration:
In this scenario we simulate 2 populations
```console
//Input parameters for fastsimcoal2
2
```
We want to calculate the effictive populaiton size 
```console
/Deme sizes (haploid number of genes)
NePOP1$
NePOP2$
```
Our sample size is the number of HAPLOID individuals, so X samples x 2. 
```console
//Sample Sizes - the HAPLOID sample size and generations since sampling (always 0)
68 0
68 0
```
We set the growth rate is set to 0
```console
//Growth rates
0
0
```
We want fastsimcoal2 to estimate the migration rate from the X to Y population and vice versa so in the migration matrix, we set the migration rate for those two scenarios as variables to use in the .est file. 
```console
//NoMig matrix 1 :
0.000 MIG0-1$ 
MIG1-0$ 0.000
```
We want to simulate divergence and migration between two populations so we add a historical event that tells the program this and asks it to estimate when the two populations diverged. 
```console
//Historical event-This estimates when the populations were isolated the identity of the source or sink does not matter here: time (in gen), source, sink, proportion of migrants, new deme size, new growth rate, new migration matrix. Below models the movement of 100% of lineages in pop0 to pop1 and estimates the number of generations since this happened
1 historical event
T_DIV$ 0 1 1 1 0 keep
```
These paramaters describe the genomic data to the program and should not be changed
```console
//Number of independent (unlinked) chromosomes, and "chromosome structure" flag:  0 for identical structure across chromosomes, and 1 for different structures on different chromosomes.
1 0
//Number of contiguous linkage blocks in chromosome 1:
1
```
Finally, we define the data types, number of loci, recomination rate, and mutation rate. The mutation rate for basidimycetes is generally not known but 1e-8 is standard.
```console
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 1e-8 OUTEXP
```

