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

# GBR - 91 genomes
gbr_merged = read.csv("~/Documents/STRs/ANALYSIS/population_research/1kg/data/GBR/merged/merged_GBR_91_genomes_1Kg_EHv3.2.2.tsv",
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")
dim(gbr_merged)
# 567  12

gbr_merged = gbr_merged %>%
  filter(gene %in% "AR") %>%
  select(gene, allele, num_samples)
gbr_merged$subpopu= rep("GBR", length(gbr_merged$gene))
dim(gbr_merged)
# 17  4

# CEU - 97 genomes
ceu_merged = read.csv("~/Documents/STRs/ANALYSIS/population_research/1kg/data/CEU/merged/merged_CEU_97_genomes_1Kg_EHv3.2.2.tsv",
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")
dim(ceu_merged)
# 577  12

ceu_merged = ceu_merged %>%
  filter(gene %in% "AR") %>%
  select(gene, allele, num_samples)
ceu_merged$subpopu= rep("CEU", length(ceu_merged$gene))
dim(ceu_merged)
# 17  4

# FIN - 86 genomes
fin_merged = read.csv("~/Documents/STRs/ANALYSIS/population_research/1kg/data/FIN/merged/merged_FIN_86_genomes_1Kg_EHv3.2.2.tsv",
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")
dim(fin_merged)
# 472  12

fin_merged = fin_merged %>%
  filter(gene %in% "AR") %>%
  select(gene, allele, num_samples)
fin_merged$subpopu= rep("FIN", length(fin_merged$gene))
dim(fin_merged)
# 12  4

eur_merged = rbind(ibs_merged,
                   tsi_merged,
                   gbr_merged,
                   ceu_merged,
                   fin_merged)
dim(eur_merged)
# 77 4