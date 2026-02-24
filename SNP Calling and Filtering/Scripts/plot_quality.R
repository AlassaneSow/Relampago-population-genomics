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
