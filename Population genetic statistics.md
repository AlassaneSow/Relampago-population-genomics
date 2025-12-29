# Calculating nucleotide diversity, fixation index, absolute divergence, and Taijama's D 
## Required files
The populations file is necessary for fst and dxy. 
## Calculating statistics
We used ```PIXY``` to calculate each statisitc as follows:     
``
pixy --stats pi fst dxy tajima_d \
--vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--window_size 10000 \
--n_cores 2
``  
# Genic diversity indicies
We found it important to calculate these statistics for only genic regions as well. In order to do so we created a bed file of all genic regions in each alignment by filtering the GTF file to only include genes  
``
``
and then retaining only the chromosome ID, start, and end and converting to bed format  
``
``
We then calculated the statistics as follows
``
pixy --stats pi fst dxy tajima_d \
--vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--bed_file genomic_windows.bed
``



