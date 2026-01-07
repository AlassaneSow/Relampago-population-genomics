# SNP calling using GATK 

## Required Files
paired.sam from bwa alignment

First need to call variants using GATK ```HaplotypeCaller```

Then we need to combine the gVCF files from HaplotypeCaller into one VCF using ```CombineGVCFs``` 

Next we need to joint call the SNPs using ```GenotypeGVCFs``` 

Use ```VariantFiltration``` default settings to filter variants

Split output by sample

Select only SNPs from output

