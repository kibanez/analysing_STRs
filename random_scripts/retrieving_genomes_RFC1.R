# Objective: retrieve the genomes (platekeys) for which EH (both v2.5.5 and v3.1.2) have estimated an insertion-STR in RFC1 or Canvas
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/CANVAS_RFC1/")

# Input merged VCF file for EHv2.5.5
canvas_v2 = read.csv("~/Documents/STRs/data/research/EH_2.5.5_research_August2019/EH_output_v2.5.5_research_August_2019/merged_genotypeUpdated/merged_loci_86457_research_genomes_new_loci_EHv2.5.5_summer2019_removingListVcfFiles.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(canvas_v2)
# 3983  11

canvas_v2 = canvas_v2 %>% filter(gene %in% "CANVAS_AAGGG")
dim(canvas_v2)
# 89  11

# We need to retrieve now the list of platekeys within `list_samples` that are separated by ';'
l_genomes_canvas_v2 = c()
for (i in 1:length(canvas_v2$list_samples)){
  aux_l_genomes = canvas_v2$list_samples[i]
  l_genomes_canvas_v2 = c(l_genomes_canvas_v2,
                          strsplit(aux_l_genomes, ';')[[1]])
  
}

l_genomes_canvas_v2 = unique(l_genomes_canvas_v2)
length(l_genomes_canvas_v2)
# 7159

write.table(l_genomes_canvas_v2,
            "./list_7159_genomes_RFC1_EHv2.5.5.tsv",
            quote = F,
            row.names = F,
            col.names = F)


# Input merged VCF file for EHv3.1.2
canvas_v3 = read.csv("~/Documents/STRs/data/research/EH_3.1.2_research_October2019/merged_genotypeUpdated/merged_loci_82565_research_genomes_EHv3.1.2_october2019.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(canvas_v3)
# 4495  12

canvas_v3 = canvas_v3 %>% filter(gene %in% "RFC1")
dim(canvas_v3)
# 248  12

# We need to retrieve now the list of platekeys within `list_samples` that are separated by ';'
l_genomes_canvas_v3 = c()
for (i in 1:length(canvas_v3$list_samples)){
  aux_l_genomes = canvas_v3$list_samples[i]
  l_genomes_canvas_v3 = c(l_genomes_canvas_v3,
                          strsplit(aux_l_genomes, ';')[[1]])
  
}

l_genomes_canvas_v3 = unique(l_genomes_canvas_v3)
length(l_genomes_canvas_v3)
# 82564

write.table(l_genomes_canvas_v3,
            "list_82564_genomes_RFC1_EHv3.1.2.tsv",
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)
