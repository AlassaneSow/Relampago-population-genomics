setwd("/blue/matthewsmith/a.sow/Sebastian_Data")
library(LEA)
library(ggplot2)

obs <- snmf("LEA_genotypes_for_strucutre.geno", K = 1:5, repetitions = 10, entropy = TRUE, project = "new", CPU = 12)

Fp_cross_entropy <- data.frame(t(summary(obs)$crossEntropy), K = 1:5, species = "Fp")

png("CrossEntropy.png", height=600, width=600)

Kplot <- ggplot(Fp_cross_entropy, aes(x=K, y=mean))+
geom_line()+
labs(x="K", y="Cross entropy")+
scale_x_continuous(breaks=1:5)

print(Kplot)
dev.off()
