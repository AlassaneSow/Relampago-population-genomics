# Relampago population genomics
[![BioProject](https://img.shields.io/badge/BioProject-PRJNAXXXXX-a03)]([https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1400107/](https://github.com/AlassaneSow/Relampago-population-genomics/edit/main/README.md))  
[![SRA](https://img.shields.io/badge/Raw_Reads-SRRXXXXX-a03)]([https://www.ncbi.nlm.nih.gov/sra/SRX31740809[accn]](https://github.com/AlassaneSow/Relampago-population-genomics/edit/main/README.md)) 
[![Mito](https://img.shields.io/badge/Mito_Genome-PXXXXXX-a03)]([https://www.ncbi.nlm.nih.gov/nuccore/PX891014.1](https://github.com/AlassaneSow/Relampago-population-genomics/edit/main/README.md))

This repository contains code and supplementary information used for the analysis of the population structure, ecology, and evolution of  _Parvodontia relampaga_. Each folder represents one major step of the analysis and in each folder you will find:
  1) A README that contains notes, supplementary figures, and summarises the commands/programs used, my thought process behind using them, and occasional trouble shooting information
  2) Scripts used for each part of the analysis

note to self. 
one paper describing population strucutre, naivte or not? host-pathogen genotyper interactions  
one paper desribing migration patterns, genes involved in local adaptation, and ecological niche modeling   
# need to make picture of flow chart bruh see https://onlinelibrary.wiley.com/doi/10.1111/mec.70321 for good flowchart
# Refer to and cite this in approriate folders https://speciationgenomics.github.io/

Note - With the availability of several hundreds to tens of thousands of loci in the population genomic realm, neutral loci can be differentiated from adaptive loci by calculating locus-statistics such as the fixation index (FST), which is a statistical measure of population differentiation based upon genetic distance. Loci with outlier FST values are likely to be under selection (adaptive loci) and hence should be removed from the analyses for making accurate inferences about demographic processes (Grünwald et al. 2016). ## Need to double check that the SNP filtering does this prior to PCA etc etc. 

# TEMPORARY CHEAT SHEET
Extract sample names from vcf like so 
```console
bcftools query -l CLEAN_ALLSITES.vcf.gz > samples.txt
```

# Demographic modeling
https://speciationgenomics.github.io/fastsimcoal2/
