library(gridExtra)
library(ggplot2)
setwd("")
SNPS <- read.table("raw_snp_scores.table", header=TRUE)

png("DP.png", height=2000, width=1500)
DP <- ggplot(SNPS, aes(x=DP))+
      geom_density(fill = "purple", alpha=0.3)+
      geom_vline(xintercept=c(10,6200))
print(DP)
dev.off()

png("QD.png", height=2000, width=1500)
QD <- ggplot(SNPS, aes(x=QD))+
      geom_density(fill = "purple", alpha=0.3)+
      geom_vline(xintercept=2, size=0.7)
print(QD)
dev.off()

png("FS.png", height=2000, width=1500)
FS <- ggplot(SNPS, aes(x=FS))+
      geom_density(fill = "purple", alpha=0.3)+
      geom_vline(xintercept=60, size=0.7)
print(FS)
dev.off()

png("MQ.png", height=2000, width=1500)
MQ <- ggplot(SNPS, aes(x=MQ))+
      geom_density(fill = "purple",alpha=0.3)+
      geom_vline(xintercept=40, size=0.7)
print(MQ)
dev.off()

png("MQRankSum.png", height=2000, width=1500)
MQRankSum <- ggplot(SNPS, aes(x=MQRankSum))+
      geom_density(fill = "purple", alpha=0.3)
print(MQRankSum)
dev.off()

png("SOR.png", height=2000, width=1500)
SOR <- ggplot(SNPS, aes(x=SOR))+
      geom_density(fill = "purple",alpha=0.3)+
      geom_vline(xintercept=4, size=0.7)
print(SOR)
dev.off()

png("snp_quality_plots.png", height=2000, width=1500)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD, DP, FS, MQ, MQRankSum, SOR, nrow=4)
dev.off()
