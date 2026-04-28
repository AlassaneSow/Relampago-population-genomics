# Variant calling using the Genome Analysis ToolKit 
This code was largely based on the methods described by [Treindl et al. 2023](https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2023.1129867/full), [Tremble et al. 2022](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.18521), and notes provided by Dmytro Kryvokhyzha in [GATK: the best practice for genotype calling in a non-model organism](https://evodify.com/gatk-in-non-model-organism/)

### Required Files
Input - BAM file containing reads   
Reference - Reference genome 
### Required Scripts
HCall.sh
## First we call variants using ```HaplotypeCaller```    

Because ```HaplotypeCaller``` can only handle one sample at a time we need to call each alignment one at a time. To do this we created a list file containing every alignment.   
```console
ls /path_to_alignments/*bam > bam_list.txt
```
We also have to index the reference genome like so
```console
samtools faidx /path/to/refernce/reference.fna
```
```console
gatk CreateSequenceDictionary \
   -R /path/to/refernce/reference.fna \
   -O /path/to/refernce/reference.dict
```

Then we ran ```HaplotypeCaller``` on all the alignments. Each variable is explained below. 
```console
REF="/path/to/reference.fna"
names="/path/to/names.txt"
OUT="path/to/snp_calling_folder"
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${names}")
BAM="/path/to/filtered_assemblies"
```
```console
gatk HaplotypeCaller \   
-R ${REF} \  
-I ${BAM} \    
-ERC GVCF \  
-O $OUT/${SAMPLE}.g.vcf.gz
```    
Then we combinened the gVCF files from HaplotypeCaller into one VCF using ```GenomicsDBImport```  
```console
gatk GenomicsDBImport \
-R /path_to_reference \
--genomicsdb-workspace-path cohort_db \
--sample-name-map gvcf_map.txt \
```
To make the gvcf_map.txt, simply run this in the HaplotypeCaller output directory  
```console
for f in *.g.vcf.gz; do 
    sample=$(basename "$f" .g.vcf.gz)
    echo -e "${sample}\t$(pwd)/${f}"
done > gvcf_map.txt
```
NOTES
```console
[a.sow@login7 jobs]$ gatk CreateSequenceDictionary -R /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna
INFO:    underlay of /etc/localtime required more than 50 (100) bind mounts
Using GATK jar /gatk/gatk-package-4.6.2.0-local.jar
Running:
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /gatk/gatk-package-4.6.2.0-local.jar CreateSequenceDictionary -R /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/blue/plp6235/asow/Fomotopsis/ref/tmp
INFO    2026-04-03 14:22:59     CreateSequenceDictionary        Output dictionary will be written in /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.dict
14:22:59.922 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/gatk/gatk-package-4.6.2.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
[Fri Apr 03 14:22:59 GMT 2026] CreateSequenceDictionary --REFERENCE /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna --TRUNCATE_NAMES_AT_WHITESPACE true --NUM_SEQUENCES 2147483647 --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 2 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false
[Fri Apr 03 14:23:00 GMT 2026] Executing as a.sow@login7.ufhpc on Linux 5.14.0-503.40.1.el9_5.x86_64 amd64; OpenJDK 64-Bit Server VM 17.0.12+7-Ubuntu-1ubuntu222.04; Deflater: Intel; Inflater: Intel; Provider GCS is available; Picard version: Version:4.6.2.0
[Fri Apr 03 14:23:00 GMT 2026] picard.sam.CreateSequenceDictionary done. Elapsed time: 0.00 minutes.
Runtime.totalMemory()=285212672
To get help, see http://broadinstitute.github.io/picard/index.html#GettingHelp
picard.PicardException: file:///blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.dict already exists.  Delete this file and try again, or specify a different output file.
        at picard.sam.CreateSequenceDictionary.doWork(CreateSequenceDictionary.java:227)
        at picard.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:281)
        at org.broadinstitute.hellbender.cmdline.PicardCommandLineProgramExecutor.instanceMain(PicardCommandLineProgramExecutor.java:37)
        at org.broadinstitute.hellbender.Main.runCommandLineProgram(Main.java:166)
        at org.broadinstitute.hellbender.Main.mainEntry(Main.java:209)
        at org.broadinstitute.hellbender.Main.main(Main.java:306)
[a.sow@login7 jobs]$ ls
calc_coverage.sh           combine_gvcf_28670586.out  combine_gvcf.sh     mapping_stats_and_filtering.sh  outputs_logs
combine_gvcf_28667660.out  combine_gvcf_28670713.out  haplotypecaller.sh  map_sam2bam.sh                  trim_reads.sh
[a.sow@login7 jobs]$ cd /blue/plp6235/asow/Fomotopsis/ref/
[a.sow@login7 ref]$ ml samtools
[a.sow@login7 ref]$ samtools faidx /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna
awk '{print $1 ":" 1 "-" $2}' /blue/plp6235/asow/Fomotopsis/ref/fomotopsis_reference.fna.fai > /blue/plp6235/asow/Fomotopsis/ref/genome_intervals.list
[a.sow@login7 ref]$ nano genome_intervals.list
```

## After calling variants, we preform joint genotyping using ```GenotypeGVCFs``` 
```console
gatk GenotypeGVCFs \
-R reference.fasta \
-V gendb://cohort_db \
--include-non-variant-sites
-O cohort.vcf.gz
```
## Next, we analyze the quality of the SNPs before filtering. 
To do this we first extract just SNPs from the cohort.vcf.gz
```console
gatk SelectVariants \
-R reference.fasta \
-V cohort.vcf.gz \
-selectType SNP \
-o cohort_all_SNPs.vcf.gz
```
Then we extract the quality metrics from the vcf.gz
```console
gatk VariantsToTable \
-R reference.fasta \
-V cohort_all_SNPs.vcf.gz \
-F CHROM -F POS -F QUAL -F QD -GF DP -F MQ -F MQRankSum -F FS -F ReadPOSRankSum -F SOR \
--allowMissingData \
-o SNP_Scores.table
```
Lastly, we plot the quality metrics in R.
```r
library('gridExtra')
library('ggplot2')

SNPs <- read.csv("path/to/SNP_Scores.table")

DP <- ggplot(SNPS, aes(x=DP, fill=purple))+
      geom_density(alpha=0.3)+
      geom_vline(xintercept=c(10,6200))

QD <- ggplot(SNPS, aes(x=QD, fill=purple))+
      geom_density(alpha=0.3)+
      geom_vline(xintercept=2, size=0.7)

FS <- ggplot(SNPS, aes(x=FS, fill=purple))+
      geom_density(alpha=0.3)
      geom_vline(xintercept=60, size=0.7)

MQ <- ggplot(SNPS, aes(x=MQ, fill=purple))+
      geom_density(alpha=0.3)+
      geom_vline(xintercept=40, size=0.7)

MQRankSum <- ggplot(SNPS, aes(x=MQRankSum, fill=purple))+
      geom_density(alpha=0.3)

SOR <- ggplot(SNPS, aes(x=SOR, fill=purple))+
      geom_density(alpha=0.3)
      geom_vline(xintercept=4, size=0.7)

ReadPosRankSum <- ggplot(SNPS, aes(x=ReadPosRankSum, fill=purple))+
      geom_density(alpha=0.3)+
      geom_vline(xintercept=c(10,-10), size=0.7)

svg("path/to/output/QualityPlots.svg", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, nrow=4)
dev.off()
```

## Now, based on the information above and GATK recommendations, we can use ```VariantFiltration```, ```vcftools```, and  ```bcftools``` to remove low quality SNPs.   
Please note that because ```pixy``` requires an invariant+variant sites file AND some of the SNP filtering commands remove invariant sites, we have to filter each site type sepratley. In practice this means we have to make SNP-only and an invariant-only files, filter them, and recombine them. We also used the SNP only filtering to make a vcf with only biallelic SNPs for our other analyses (e.g., ADMIXTURE).

### General filtering
```console
vcftools --gzvcf my_vcf.vcf.gz \
--remove-indels \
--max-missing 0.95 \
--min-meanDP 10 \
--max-meanDP 75 \
--recode --stdout | bgzip -c > missingness_filtered_ALLSITES.vcf.gz
```
### Filtering **SNPs ONLY**
We extract just the SNPs using vcftools
```console
vcftools --gzvcf missingness_filtered_ALLSITES.vcf.gz \
--mac 1 \
--recode --stdout | bgzip -c > missingness_filtered_SNPS.vcf.gz
```
And then filter them using ```VariantFiltration```
```console
gatk VariantFiltration \
-R ${REF} \
-V missingness_filtered_SNPS.vcf.gz \
-O quality_labeled_SNPS.vcf.gz \
--filter "MQ < 40.0" --filter-name "MQ" \
--filter "FS > 60.0" --filter-name "FS" \
--filter "QD < 2" --filter-name "QD"

bcftools view -f PASS quality_labeled_SNPS.vcf.gz -Oz -o quality_filtered_SNPS.vcf.gz
```
This data can now be combined with the invariant sites file shown below!  
Most population genetic analyses require bi-allelic SNPs, to retain only bi-allelic SNPs we ran
```
vcftools --gzvcf quality_filtered_SNPS.vcf.gz \
--min-alleles 2 \
--max-alleles 2 \
--recode --stdout | bgzip -c > quality_filtered_BIALLELIC_SNPS.vcf.gz
```
Optional removal of non-chromosomal sites. Run
```console
bcftools view \
-r Chr1,Chr2..... \
filtered_SNPS.vcf.gz -Oz -o nuclear_filtered_SNPS.vcf.gz #In the command above you can use -R and give a list of chromosome names rather than typing them out.
```
We started with X SNPs
```console
grep -vc "^#" missingness_filtered_SNPS.vcf.gz
```
X SNPs are included in our PIXY analysis
```console
grep -vc "^#" quality_filtered_SNPS.vcf.gz
```
and there are X high-quality bi-allelic SNPs
```console
grep -vc "^#" quality_filtered_BIALLELIC_SNPS.vcf.gz
```

### Filtering **invariant sites**
Because 1) the quality filters used in the SNP section do not apply to invariant sites and 2) we already removed sites with missing data or low/high depth we simply extract invariant sites from the vcf and combine with the filtered SNPs.
```console
gatk SelectVariants \
  -R ${REF} \
  -V missingness_filtered_ALLSITES.vcf.gz \
  --select-type-to-include NO_VARIATION \
  -O quality_filtered_INVARIANT.vcf.gz
```
### Combine, index, and sort files
```console
tabix -p vcf quality_filtered_INVARIANT.vcf.gz
tabix -p vcf quality_filtered_SNPS.vcf.gz

bcftools concat \
--allow-overlaps \
quality_filtered_INVARIANT.vcf.gz quality_filtered_SNPS.vcf.gz \
-O z -o quality_filtered_ALLSITES.vcf.gz

bcftools sort quality_filtered_ALLSITES.vcf.gz -Oz -o pixy_input_quality_filtered_ALLSITES.vcf.gz
bcftools index pixy_input_quality_filtered_ALLSITES.vcf.gz
```
## Now we have a file containing bi-allelic SNPs that we can use for ADMIXTURE and PCA and a file with invariant+variant sites for PIXY. With this data we first [linkage-disequillibrium pruned](/Linkage%20Disequilibrium/README.md) the SNPs before analyzing the [population strucutre](/Population%20Structure/Population%20Structure.md) of _Parvodontia relampaga_. 


