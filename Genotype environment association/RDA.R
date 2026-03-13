#RDA
library(adegenet)
library(vcfR)
library(tidyverse)
library(caret)
library(vegan)
#Because RDA can not handle missing data, we need to perform mean imputation of missing alleles.
#First, we create a matrix so we can modify the genlight object then we re-code the values for each allele
geno_matrix <- as.matrix(data)  # individuals x SNPs
geno_matrix[geno_matrix == 0] <- -1  # recode homozygote ref
geno_matrix[geno_matrix == 1] <- 0   # recode heterozygote
geno_matrix[geno_matrix == 2] <- 1   # recode homozygote alt
#Then, we use mean imputation to make NAs the mean genotype value at each locus.   
geno_matrix[is.na(geno_matrix)] <- apply(geno_matrix, 2, function(x) mean(x, na.rm = TRUE))[col(geno_matrix)[is.na(geno_matrix)]]
#Lastly we assign new row names to the matrix. The syntax for gsub() is gsub(pattern, replacement,text to modify). Here, | means OR. The code below removes file paths and .bam from the names              
rownames(geno_matrix) <- gsub("/media/hdd2tb/popgen/Pilze/PopGen/BAM-Ac/|\\.bam|'|/media/hdd2tb/popgen/Pilze/PopGen/BAM-Fp/|\\.bam|","",rownames(geno_matrix))
#Next, load a data frame with metadata and environmental variables. 
metadata <- readRDS("env.rds")
#Here, the bioclimatic variables are extracted from each coordinate and merged with the rest of the metadata. 
#For my data I need to first extract the bio-climatic variables from the worldclim data
#CODE HERE
#Then merge with my metadata file. 
#CODE HERE. 
#Because this already done by Sebastian et al. we can just proceed. 
#We reorder the rows  to match those in the geno_matrix and remove those that don't have same Samples
metadata <- metadata %>%
  slice(match(rownames(geno_matrix), Sample))
#double check that the row names and sample IDs are the same
rownames(geno_matrix) == metadata$sample_name
#optional removal of whitspaces
metadata$count <- gsub("\\s", "", metadata$count)
#Next we extract the bioclim variables from the metadata dataframe, and remove the variables that are correlated >0.7
bioclim <- metadata[,25:43]
correlations <- cor(bioclim)
High_correlations <- findCorrelation(correlations, cutoff=0.7, exact = TRUE)
bioclim[,High_correlations] <- list(NULL)
#Running the RDA will take a while so we should run and analyze it on the Hipergator. Save the genotype matrix and climate variables
saveRDS(bioclim,"bioclim_variables_for_rda.rds")
saveRDS(geno_matrix, "genotype_matrix_for_rda.rds")
rda_model <- rda(geno_matrix ~., data=bioclim)
saveRDS(rda_model, "rda_model.rds")
RsquareAdj(rda_model)
#About 5% of the variance was explained by the bioclim variables. This is expected because we would assume that most SNPs are neutral mutations
#Now can run an anova test for the whole model. The anova.cca() command is made for RDA, CCA, or drBDA.
rda_model_fit <- anova.cca(rda_model)
#Using by="terms" will test the signifcance of each term/variable and by="axis" will test the significane of each constrained axis.
rda_variable_fit <- anova.cca(rda_model, by = "terms", permutations = how(nperm = 99),
                        parallel = 10)

#       Df Variance      F Pr(>F)
#bio_1     1     9385 1.1361   0.04 *
#bio_2     1     9301 1.1259   0.06 .
#bio_8     1    10768 1.3034   0.01 **
#bio_9     1     9539 1.1547   0.08 .
#bio_10    1    16010 1.9380   0.01 **
#bio_15    1     8506 1.0296   0.23
#bio_18    1     8406 1.0176   0.29
#Residual 31   256094
#The anova test tells us that bio_8,bio_10,and bio_1 are significanlty associated with SNP variation.
# Next, we extract eigenvalues from the constrained and unconstrained axes
eig_values <- rda_model$CCA$eig  # constrained axes
eig_total  <- rda_model$CA$eig   # unconstrained (residual) axes
prop_explained <- eig_vals / sum(rda_model$CCA$eig)
#The proportion of variance explained by each axis is as follows
#RDA1      RDA2      RDA3      RDA4      RDA5      RDA6      RDA7
#0.3152886 0.1236675 0.1223969 0.1197129 0.1113647 0.1042588 0.1033106

#The authors did not scale the bioclim variables. What happens if I do? 

#Plot individuals
#Extract information to make the bioclim vectors
envrionmental_vectors <- scores(rda_model, display = "bp", scaling = 2) #extract coordinates
scaled_environmental_vectors <- envrionmental_vectors*12 #scale for readability
df_envrionmental_vectors <- as.data.frame(scaled_environmental_vectors) 
df_envrionmental_vectors$Variable <- rownames(scaled_environmental_vectors) #add variable column
#Extract site coordinates for each individual 
sp <- scores(rda_model, display = "sites", scaling = 2)
df_sp <- as.data.frame(sp)
#Plot!  
ggplot()+ 
  geom_point(data=df_sp, aes(x=RDA1, y=RDA2), size=3)+ #specimen points
  geom_segment(data=df_environmental_vectors,aes(x=0,y=0,xend=RDA1,yend=RDA2), arrow=arrow(length=unit(0.25, "cm")))+ #plot environmental vectors
  geom_text(data=df_envrionmental_vectors, aes(x=RDA1, y=RDA2, label=Variable), size=4)+ #labels for arrows
  labs(x="RDA1", y="RDA2")+
  theme_minimal()
# Note to self to also figure out how to color by population
#Plot SNPs
