library(tidyverse)

fst <- read_tsv("path/to/pixy_fst_output.tsv")
dxy <- read_tsv("path/to/pixy_dxy_output.tsv")
genes <- read_tsv("path/to/genes.bed",
                  col_names = c("chromosome","start","end","gene_id"))
data <- fst %>%
  select(chromosome,start,end,pop1,pop2,avg_wc_fst) %>%
  left_join(
  dxy %>%
    select(chromsome, start,end,pop1,pop2,avg_dxy),
  by = c("chromsome", "start","end","pop1","pop2","avg_dxy")
)
data <- data %>%
  left_join(genes,
            by = c("chromosome","start","end"))
divergent_loci <- data %>%
  group_by(pop1,pop2) %>%
  filter(
    avg_wc_fst >= quantile(avg_wc_fst,0.95,na.rm=TRUE),
    avg_dxy >= quantile(avg_dxy, 0.95, na.rm=TRUE)) %>%
    ungroup()
write_tsv(data, "raw_divergent_loci.tsv")
write_tsv(divergent_loci, "divergent_loci_fst_dxy.tsv")
