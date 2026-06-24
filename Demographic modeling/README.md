# Demographic modeling using fastsimcoal2
This code is based on the [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/fastsimcoal2/), the [fastsimcoal2 manual](https://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal28.pdf), and methods used by Tremble et al. [2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521).
## Preparing data
Before running ```fastsimcoal2``` we have to estimate a site-frequency spectrum using easy SFS. To do this we will need the LD-pruned vcf and a populations file. We will use the popualtions identified in the K=X [ADMIXTURE plots](https://github.com/AlassaneSow/Relampago-population-genomics/tree/main/Population%20Structure#structure).  
  
Although the vcf contains very high quality SNPs, ```easysfs``` can not use missing data and we may have to downsample our data to retain the highest number of sites. If there is lots of missing data or low covergae, the alternate program ANGS can create an SDS and account for genotype uncertanties.
```console
# estimate the best projections
easSFS -i ${VCF} -p ${POP} -a -f --preview
# then create SFS files
```

