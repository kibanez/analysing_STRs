# Objective: analyse 198 loci annalysed by Zhongbo at UCL, specific batch december 2020
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
library(purrr); packageDescription ("purrr", fields = "Version") #"0.3.3"

# Set working dir
setwd("~/Documents/STRs/data/research/batch_december2020/output_EHv3.2.2_vcfs/merged/")

# Functions
source("~/git/analysing_STRs/functions/plot_together_histo_boxplot.R")
source("~/git/analysing_STRs/functions/computing_percentiles.R")

# load merged august data
merged_data = read.csv("./merged_93430_genomes_EHv322_batch_december2020.tsv",
                       stringsAsFactors = F,
                       sep = "\t",
                       header = T)
dim(merged_data)
# 6612  12

# Some genes have `/` and R can interpret as part of a path within the system
merged_data$gene = gsub("/", "_", merged_data$gene)

# load clinical data - changing to RE V11 (since we are sharing with external groups)
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_301220_V11.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 2119961  36

# Let's keep only germline genomes
clin_data = clin_data %>%
  filter(grepl("germline", type))
dim(clin_data)
# 2092419 36

# Checking which genomes we have
table(clin_data$type)
#cancer germline experimental germline rare disease germline 
#24649                   84               2067686

# List of platekeys corresponding to ONLY PROBANDS in Neuro
df_only_probands = clin_data %>%
  filter(is.na(biological_relationship_to_proband) |
           biological_relationship_to_proband %in% "N/A" | 
           biological_relationship_to_proband %in% "Proband" |
           programme %in% "Cancer")

l_platekeys_probands_neuro = df_only_probands %>%
  filter(grepl("neuro", disease_group, ignore.case = TRUE)) %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_platekeys_probands_neuro)
# 13796

# There are some platekeys (16k) that have ',', which means that PID is associated with more than one platekey
l_platekeys_probands_neuro_unique = c()
for (i in 1:length(l_platekeys_probands_neuro)){
  if (grepl(',',l_platekeys_probands_neuro[i])){
    list_platekeys = strsplit(l_platekeys_probands_neuro[i], ",")[[1]]
    list_platekeys = gsub(" ", "", list_platekeys, fixed = TRUE)
    l_platekeys_probands_neuro_unique = c(l_platekeys_probands_neuro_unique,
                                    max(list_platekeys))
  }else{
    l_platekeys_probands_neuro_unique = c(l_platekeys_probands_neuro_unique,
                                    l_platekeys_probands_neuro[i])
  }
}
length(l_platekeys_probands_neuro_unique)
# 13796

# List of platekeys corresponding to ONLY PROBANDS but NOT in Neuro
# First probands
df_only_probands_notNeuro = df_only_probands %>%
  filter(!grepl("neuro", disease_group, ignore.case = TRUE))
dim(df_only_probands_notNeuro)
# 687859 36

l_platekeys_probands_notNeuro = df_only_probands_notNeuro %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_platekeys_probands_notNeuro)
# 35476

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
# 35476

# 1. Merge GRCh37 and GRCh38 info, since chromosome names are different
# GRCh38 are chr1, chr2, chr3 while GRCh37 are 1,2,3
# In SHARP everything should be GRCh38, because Andy sent us coordinates only in GRCh38
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

# Let's focus on SHARP genes
l_genes = unique(merged_data$gene)
length(l_genes)
# 198

# Output folder
output_folder = 'EHv322_batch_december2020'
#dir.create(output_folder)

for(i in 1:length(l_genes)){
  plot_together_histo_boxplot(df_input = merged_data,
                              gene_name = l_genes[i],
                              output_folder = output_folder,
                              l_platekeys_probands_neuro_unique,
                              l_platekeys_probands_notNeuro_unique)
}

# Summarise report with quantiles for all genes
# Based on 93k genomes
df_percentiles = computing_percentiles(merged_data)

# Write down first a line with total number of probands in neuro and probands in not neuro
to_write = cbind("93430",
                 length(l_platekeys_probands_neuro_unique),
                 length(l_platekeys_probands_notNeuro_unique))
to_write = rbind(c("Total number of genomes","Total number of probands recruited under Neurological disorders", "Total number of probands NOT recruited under Neurological disorders"),
                 to_write,
                 cbind("", "", ""))
write.table(to_write, 
            "./EHv322_batch_december2020/summary_stats_quantiles_198_Zhongbo_genes.tsv", 
            sep='\t', 
            quote = F,
            row.names=F, 
            col.names=F)
write.table(df_percentiles,
            "./EHv322_batch_december2020/summary_stats_quantiles_198_Zhongbo_genes.tsv", 
            append = T, 
            sep = "\t", 
            quote = F,
            row.names = F,
            col.names = T)