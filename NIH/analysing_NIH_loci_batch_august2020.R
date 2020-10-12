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
library(magick); packageDescription ("magick", fields = "Version") #"2.4.0"
library(cowplot); packageDescription ("cowplot", fields = "Version") #"1.0.0"
library(ggpubr); packageDescription ("ggpubr", fields = "Version") #"1.0.0"


# Set working dir
setwd("~/Documents/STRs/ANALYSIS/NIH/")

# Functions
source("~/git/analysing_STRs/functions/plot_together_histo_boxplot.R")

# load merged august data
merged_data = read.csv("~/Documents/STRs/data/research/batch_august2020/output_EHv3.2.2_vcfs/merged/merged_93446_genomes_EHv322_batch_august2020.tsv",
                       stringsAsFactors = F,
                       sep = "\t",
                       header = T)
dim(merged_data)
# 27238  12

# load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_research_cohort/clinical_data_research_cohort_93614_PIDs_merging_RE_V1toV10.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 152337  32

# List of platekeys corresponding to ONLY PROBANDS
df_only_probands = clin_data %>%
  filter(is.na(biological_relationship_to_proband) |
           biological_relationship_to_proband %in% "N/A" | 
           biological_relationship_to_proband %in% "Proband" |
           programme %in% "Cancer")

l_platekeys_probands = df_only_probands %>%
  select(list_platekeys1) %>%
  unique() %>%
  pull()
length(l_platekeys_probands)
# 52191

# There are some platekeys (16k) that have ',', which means that PID is associated with more than one platekey
l_platekeys_probands_unique = c()
for (i in 1:length(l_platekeys_probands)){
  if (grepl(',',l_platekeys_probands[i])){
    list_platekeys = strsplit(l_platekeys_probands[i], ",")[[1]]
    list_platekeys = gsub(" ", "", list_platekeys, fixed = TRUE)
    l_platekeys_probands_unique = c(l_platekeys_probands_unique,
                                    max(list_platekeys))
  }else{
    l_platekeys_probands_unique = c(l_platekeys_probands_unique,
                                    l_platekeys_probands[i])
  }
}
length(l_platekeys_probands_unique)
# 52191

# List of platekeys corresponding to ONLY PROBANDS but NOT in Neuro
# First probands
df_only_probands_notNeuro = df_only_probands %>%
  filter(!grepl("neuro", list_disease_group, ignore.case = TRUE))
dim(df_only_probands_notNeuro)
# 63266 32

l_platekeys_probands_notNeuro = df_only_probands_notNeuro %>%
  select(list_platekeys1) %>%
  unique() %>%
  pull()
length(l_platekeys_probands_notNeuro)
# 37701

# There are some platekeys (16k) that have ',', which means that PID is associated with more than one platekey
l_platekeys_probands_notNeuro_unique = c()
for (i in 1:length(l_platekeys_probands_notNeuro)){
  if (grepl(',',l_platekeys_probands_notNeuro[i])){
    list_platekeys = strsplit(l_platekeys_probands_notNeuro[i], ",")[[1]]
    list_platekeys = gsub(" ", "", list_platekeys, fixed = TRUE)
    l_platekeys_probands_notNeuro_unique = c(l_platekeys_probands_notNeuro_unique,
                                             max(list_platekeys))
  }else{
    l_platekeys_probands_notNeuro_unique = c(l_platekeys_probands_notNeuro_unique,
                                             l_platekeys_probands_notNeuro[i])
  }
}
length(l_platekeys_probands_notNeuro_unique)
# 37701

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

l_nih = unique(l_genes[which(grepl("^chr", l_genes, ignore.case = TRUE))])
length(l_nih)
# 192

nih_merged_data = merged_data %>%
  filter(gene %in% l_nih)
dim(nih_merged_data)
# 12608  13

# Output folder
output_folder = 'EHv322_batch_august2020'
#dir.create(output_folder)

for(i in 1:length(l_nih)){
  plot_together_histo_boxplot(df_input = nih_merged_data,
                              gene_name = l_nih[i],
                              output_folder = output_folder)
  
}  