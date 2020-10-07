# Objective: analyse loci analysed by NIH within 100kGP, after EHdn 
# For EHv322
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.1"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/NIH/")

# Function
source("~/git/analysing_STRs/functions/plot_gene_mergingAssemblies.R")

# load merged august data
merged_data = read.csv("~/Documents/STRs/data/research/batch_august2020/output_EHv3.2.2_vcfs/merged/merged_93446_genomes_EHv322_batch_august2020.tsv",
                       stringsAsFactors = F,
                       sep = "\t",
                       header = T)
dim(merged_data)
# 27238  12

# 1. Merge GRCh37 and GRCh38 info, since chromosome names are different
# GRCh38 are chr1, chr2, chr3 while GRCh37 are 1,2,3
merged_data$chr = recode(merged_data$chr,
                "1" = "chr1",
                "2" = "chr2",
                "3" = "chr3",
                "4" = "chr4",
                "5" = "chr5",
                "6" = "chr6",
                "7" = "chr7",
                "8" = "chr8",
                "9" = "chr9",
                "10" = "chr10",
                "11" = "chr11",
                "12" = "chr12",
                "13" = "chr13",
                "14" = "chr14",
                "15" = "chr15",
                "16" = "chr16",
                "17" = "chr17",
                "18" = "chr18",
                "19" = "chr19",
                "20" = "chr20",
                "21" = "chr21",
                "22" = "chr22",
                "X" = "chrX")
merged_data = merged_data %>%
  group_by(chr, gene, allele) %>%
  mutate(total_num_samples = sum(num_samples)) %>%
  ungroup() %>%
  as.data.frame() 

merged_data_simpl = merged_data %>% 
  select(chr, gene, allele, total_num_samples)
merged_data_simpl = unique(merged_data_simpl)
dim(merged_data_simpl)
# 21013  4

# Let's focus on NIH loci - the ones starting by `^chr`
l_genes = unique(merged_data$gene)
length(l_genes)
# 329

l_nih = l_genes[which(grepl("^chr", l_genes, ignore.case = TRUE))]
length(l_nih)
# 192

# Output folder
output_folder = 'EHv322_batch_august2020'
dir.create(output_folder)

for (i in 1:length(l_nih)){
  plot_gene_mergingAssemblies(merged_data_simpl, l_nih[i], output_folder)
}

