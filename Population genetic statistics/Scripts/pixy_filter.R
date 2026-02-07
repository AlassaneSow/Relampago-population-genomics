data <- read.table("path/to/pixy_output.tsv", header=TRUE, sep="\t")
fst_cut <- quantile(data$fst, 0.95, na.rm=TRUE)
dxy_cut <- quantile(data$dxy, 0.95, na.rm=TRUE)
out <- d[data$fst >= fst_cut & data$dxy >= dxy_cut, ]
write.table(out, "path/to/pixy/output/5perc_genic.tsv, sep="\t", quote=FALSE, row.names=FALSE)
