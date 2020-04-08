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
setwd("~/Documents/STRs/ANALYSIS/AR_kennedy/1Kg/")

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

sum_total = sum(ibs_merged$num_samples)
ibs_merged = ibs_merged %>%
  group_by(allele) %>%
  mutate(percent_subpopu = 100*(num_samples/sum_total)) %>%
  ungroup() %>%
  as.data.frame()
dim(ibs_merged)
# 14  5

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

sum_total = sum(tsi_merged$num_samples)
tsi_merged = tsi_merged %>%
  group_by(allele) %>%
  mutate(percent_subpopu = 100*(num_samples/sum_total)) %>%
  ungroup() %>%
  as.data.frame()
dim(tsi_merged)
# 17  5

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

sum_total = sum(gbr_merged$num_samples)
gbr_merged = gbr_merged %>%
  group_by(allele) %>%
  mutate(percent_subpopu = 100*(num_samples/sum_total)) %>%
  ungroup() %>%
  as.data.frame()
dim(gbr_merged)
# 17  5

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

sum_total = sum(ceu_merged$num_samples)
ceu_merged = ceu_merged %>%
  group_by(allele) %>%
  mutate(percent_subpopu = 100*(num_samples/sum_total)) %>%
  ungroup() %>%
  as.data.frame()
dim(ceu_merged)
# 17  5

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

sum_total = sum(ceu_merged$num_samples)
fin_merged = fin_merged %>%
  group_by(allele) %>%
  mutate(percent_subpopu = 100*(num_samples/sum_total)) %>%
  ungroup() %>%
  as.data.frame()
dim(fin_merged)
# 12  5

eur_merged = rbind(ibs_merged,
                   tsi_merged,
                   gbr_merged,
                   ceu_merged,
                   fin_merged)
dim(eur_merged)
# 77 5

# Plot percentage of each allele within each subpopulation
violin_plot = ggplot(eur_merged, aes(x = subpopu, y=allele, fill = subpopu)) +
  geom_violin() +
  coord_flip() +
  xlab("Sub-population across EUR in 1Kg dataset") + 
  ylab("Repeat sizes (repeat units)") + 
  ggtitle("AR across EUR sub-populations within 1Kg") +
  geom_boxplot(width=0.1) +
  guides(fill=guide_legend(reverse=TRUE)) 

violin_output_png = "./plots/violin_plots_across_EUR_subpopus.png"
png(violin_output_png, units="in", width=5, height=5, res=300)
print(violin_plot)
dev.off()

# Plot joint ancestry repeat-size distributions
min_value = min(eur_merged$allele, 34)
max_value = max(eur_merged$allele, 38)
joint_plot = ggplot(eur_merged, aes(x = allele, y = percent_subpopu, group = subpopu, color = subpopu)) + 
  geom_line() + 
  ylab("Allele frequency") + 
  xlab("Repeat sizes (repeat units)") + 
  geom_vline(xintercept = 34, colour = 'blue', lty = 2) + 
  geom_vline(xintercept = 38, colour = 'red', lty = 2) + 
  coord_cartesian(xlim = c(min_value,max_value))

joint_ancestries = "./plots/joint_subpopus_in_AR_1Kg.png"
png(joint_ancestries, units="in", width=5, height=5, res=300)
print(joint_plot)
dev.off()
