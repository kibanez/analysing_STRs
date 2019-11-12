# Objective: take advantage from the work GEL bioinfo research group has done estimating ancestry for ~59,356 genomes
# and analyse and study the STR distribution across them
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"

setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/")

popu_table_enriched = read.csv("./population_info_enriched_59356_by_031019.tsv",
                               header = T,
                               sep = "\t",
                               stringsAsFactors = F)
dim(popu_table_enriched)
# 59356  20

# Functions
# Function that plots the STR repeat-size frequencies for a gene/locus across the cohort
plot_gene <- function(df_input, gene_name, gene_data_normal, gene_data_pathogenic, output_folder, assembly) {
  
  threshold_normal = gene_data_normal %>% filter(grepl(gene_name, locus)) %>% select(threshold) %>% unlist() %>% unname()
  threshold_pathogenic = gene_data_pathogenic %>% filter(grepl(gene_name, locus)) %>% select(threshold) %>% unlist() %>% unname()
  
  df_gene = df_input %>% filter(gene %in% gene_name)
  
  alt_number = df_gene$allele
  
  #df_gene_barplot = data.frame(number_repeats = alt_number, af = df_gene$AF)
  df_gene_barplot = data.frame(number_repeats = alt_number, af = df_gene$num_samples)
  
  # order by 'number of repetition'
  df_gene_barplot = unique(df_gene_barplot[order(df_gene_barplot[,1]),])
  
  rownames(df_gene_barplot) = df_gene_barplot$number_repeats
  
  df_gene_barplot$number_repeats = as.numeric(df_gene_barplot$number_repeats)
  
  gene_name = paste(gene_name, assembly, sep = '_')
  pdf_name = paste(output_folder, gene_name, sep = "/")
  pdf_name = paste(pdf_name, 'pdf', sep = ".")
  
  min_value = min(df_gene_barplot$number_repeats)
  max_value = max(threshold_pathogenic + 1, df_gene_barplot$number_repeats)
  
  aux_plot = ggplot(unique(df_gene_barplot), aes(x = number_repeats, y = af)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label=af), vjust=0, size = 4, colour = "grey") +
    ylab("Allele frequency") + 
    xlab("Repeat sizes (repeat units)") + 
    ggtitle(gene_name) + 
    geom_vline(xintercept = threshold_normal, colour = 'blue', lty = 2) + 
    geom_vline(xintercept = threshold_pathogenic, colour = 'red', lty = 2) + 
    coord_cartesian(xlim = c(min_value,max_value))
  
  
  pdf(pdf_name)
  print(aux_plot)	
  dev.off()
  
}


# STR annotation, threshold including the largest normal and the smallest pathogenic sizes reported
gene_annotation_normal = '/Users/kibanez/git/analysing_STRs/threshold_largest_normal_reported_research.txt'
gene_data_normal = read.table(gene_annotation_normal, stringsAsFactors=F, header = T)

gene_annotation_pathogenic = '/Users/kibanez/git/analysing_STRs/threshold_smallest_pathogenic_reported_research.txt'
gene_data_pathogenic = read.table(gene_annotation_pathogenic, stringsAsFactors=F, header = T)



ggplot(data=dat, aes(x=peddy_pc2, y=peddy_pc1)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc2, y=peddy_pc1, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc3, y=peddy_pc1, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc4, y=peddy_pc1, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc3, y=peddy_pc2, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc4, y=peddy_pc2, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc3, y=peddy_pc4, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=1, y=peddy_ancestry_prob)) +
  geom_boxplot()

ggplot(data=dat, aes(x=peddy_ancestry_pred, y=peddy_ancestry_prob)) +
  geom_boxplot()
