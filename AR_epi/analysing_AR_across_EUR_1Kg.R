# Objective: analyse AR epi across EUR genomes in the 1Kg project
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/AR_kennedy/")

# Load data
# IBS - 102 genomes
ibs_merged = read.csv("~/Documents/STRs/ANALYSIS/population_research/1kg/data/IBS/merged/merged_IBS_102_genomes_1Kg_EHv3.2.2.tsv",
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")
dim(ibs_merged)
# 610  12

ibs_merged = ibs_merged %>%
  filter(gene %in% "AR") %>%
  select(gene, allele, num_samples)
ibs_merged$subpopu= rep("IBS", length(ibs_merged$gene))
dim(ibs_merged)
# 14  4

# TSI - 104 genomes
tsi_merged = read.csv("~/Documents/STRs/ANALYSIS/population_research/1kg/data/TSI/merged/merged_TSI_104_genomes_1Kg_EHv3.2.2.tsv",
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")
dim(tsi_merged)
# 648  12

tsi_merged = tsi_merged %>%
  filter(gene %in% "AR") %>%
  select(gene, allele, num_samples)
tsi_merged$subpopu= rep("TSI", length(tsi_merged$gene))
dim(tsi_merged)
# 17  4

