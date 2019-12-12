# Objective: take advantage from the work GEL bioinfo research group has done estimating ancestry for ~59,356 genomes
# and analyse and study the forensic STR distribution across them
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"
library(scatterplot3d); packageDescription("scatterplot3d", fields = "Version") # 0.3-41
library(ggpubr); packageDescription("ggpubr", fields = "Version") # 0.2.3


# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/HipSTR/output_HipSTR/HipSTR_output_58971_vcfs/")

# Functions
source("/Users/kibanez/git/analysing_STRs/functions/plot_violin_ancestry_without_thresholds.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene_joint_ancestries_without_thresholds.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene_without_thresholds.R")

# Create a specific folder for the figures
output_folder = "./figures_forensics/"

# Load HipSTR STR output data - merged TSV files
df_asi = read.csv('./ASI/merged/merged_forensics_loci_5947_ASI_HipSTRv0.6.2.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T)
dim(df_asi)
# 864  9

df_eas = read.csv('./EAS/merged/merged_forensics_loci_400_EAS_HipSTRv0.6.2.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_eas)
# 403  9

df_afr = read.csv('AFR/merged/merged_forensics_loci_1777_AFR_HipSTRv0.6.2.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_afr)
# 666  9

df_amr = read.csv('AMR/merged/merged_forensics_loci_797_AMR_HipSTRv0.6.2.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_amr)
# 559  9

df_eur = read.csv('EUR/merged/merged_forensics_loci_46883_EUR_HipSTRv0.6.2.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_eur)
# 1750   9


df_asi = df_asi %>% mutate(population = "ASI")
df_eas = df_eas %>% mutate(population = "EAS")
df_afr = df_afr %>% mutate(population = "AFR")
df_amr = df_amr %>% mutate(population = "AMR")
df_eur = df_eur %>% mutate(population = "EUR")

df_all = rbind(df_afr,
               df_amr,
               df_eur,
               df_eas,
               df_asi)
dim(df_all)
# 4242  10   


l_loci = sort(unique(df_all$gene))
for (i in 1:length(l_loci)){
  # Each locus - Individually
  plot_gene_without_cutoff(df_afr, l_loci[i], output_folder, "GRCh38", "AFR")
  plot_gene_without_cutoff(df_amr, l_loci[i], output_folder, "GRCh38", "AMR")
  plot_gene_without_cutoff(df_eur, l_loci[i], output_folder, "GRCh38", "EUR")
  plot_gene_without_cutoff(df_eas, l_loci[i], output_folder, "GRCh38", "EAS")
  plot_gene_without_cutoff(df_asi, l_loci[i], output_folder, "GRCh38", "ASI")
  
  # Jointly - distribution
  plot_gene_joint_ancestries_without_cutoff(df_all, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
  
  # Jointly - Violing plots
  plot_violin_ancestry_without_cutoff(df_all, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
}