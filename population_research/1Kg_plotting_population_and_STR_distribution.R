# Objective: analyse the repeat-sizes across sub-population and super-population within 1Kg cohort
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
library(tidyverse)


# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/1Kg/")

# Load 1Kg population index data
popu_info = read.csv("./integrated_call_samples_v2.20130502.ALL.ped",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)
dim(popu_info)
# 3691  17

# Load thresholds
# STR annotation, threshold including the largest normal and the smallest pathogenic sizes reported
gene_annotation_normal = '/Users/kibanez/git/analysing_STRs/threshold_largest_normal_reported_research.txt'
gene_data_normal = read.table(gene_annotation_normal, stringsAsFactors=F, header = T)

gene_annotation_pathogenic = '/Users/kibanez/git/analysing_STRs/threshold_smallest_pathogenic_reported_research.txt'
gene_data_pathogenic = read.table(gene_annotation_pathogenic, stringsAsFactors=F, header = T)


# Functions
source("/Users/kibanez/git/analysing_STRs/functions/plot_violin_ancestry_1Kg.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_violin_ancestry.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene_joint_ancestries_1Kg.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene_joint_ancestries.R")

# Load EHv3.2.2 STR merged data for each sub-population
df_merged = data.frame()
l_popus = unique(popu_info$Population)

# Remove CHD from `l_popus`
l_popus = l_popus[-19]

# Define super-populations
l_superpopu = c("AFR", "AMR", "EAS", "EUR", "SAS")

# Define sub-population and super-population
superpopulations = c("AFR","AFR","AFR","AFR","AFR","AFR","AFR", 
                     "AMR", "AMR","AMR","AMR",
                     "EUR","EUR","EUR","EUR","EUR", 
                     "EAS","EAS","EAS","EAS","EAS",
                     "SAS","SAS","SAS","SAS","SAS")
sub_populations = c("ESN", "YRI", "GWD", "LWK", "MSL", "ACB", "ASW", 
                    "MXL", "PUR", "PEL", "CLM",
                    "IBS", "TSI", "GBR", "CEU", "FIN",
                    "JPT", "CHS", "CHB", "CDX", "KHV",
                    "PJL", "STU", "BEB", "ITU", "GIH")

popu_1kg = data.frame(cbind(superpopulations, sub_populations))
popu_1kg$superpopulations = as.character(popu_1kg$superpopulations)
popu_1kg$sub_populations = as.character(popu_1kg$sub_populations)

for (i in 1:length(l_popus)){
  popu_aux = paste("~/Documents/STRs/ANALYSIS/population_research/1Kg/data/", l_popus[i] ,sep = "")
  file_aux = list.files(paste(popu_aux, "merged", sep = "/"))
  file_aux = paste(paste(popu_aux, "merged", sep = "/"), file_aux, sep = "/")
  
  df_aux = read.csv(file_aux,
                    sep  = "\t",
                    stringsAsFactors = F,
                    header = T)
  
  index_subpopu = match(l_popus[i], popu_1kg$sub_populations)
  superpopu = popu_1kg$superpopulations[index_subpopu]
  
  df_aux = df_aux %>% 
    mutate(population = l_popus[i],
           superpopulation = superpopu)
  
  df_merged = rbind(df_merged,
                    df_aux)
  
}

dim(df_merged)
# 14107  14

# Population enriched genomes are only GRCh38, we will ignore then GRCh37
output_folder = "./figures/"

l_loci = sort(unique(df_merged$gene))
# Let's focus first on the important 4 loci
l_loci = c("AR", "ATN1", "HTT", "FXN")
for (i in 1:length(l_loci)){

  # Specifying sub-population  
  for (j in 1:length(l_popus)){
    # Each locus - Individually
    plot_gene(df_merged %>% filter(population %in% l_popus[j]), l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", l_popus[j])  
  }
  
  # Jointly - distribution
  plot_gene_joint_ancestries_1Kg(df_merged, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
  
  # Jointly - Violin plots
  plot_violin_ancestry_1Kg(df_merged, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
}


# Plotting here joint ancestries per super-population (repeat-sizes across sub-populations)
l_superpopus = unique(superpopulations)
for (i in 1:length(l_loci)){
  
 # Jointly - distribution PER EACH SUPER POPULATION
  for (j in 1:length(l_superpopus)){
    plot_gene_joint_ancestries(df_merged %>% filter(superpopulation %in% l_superpopus[j]), 
                               l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, l_superpopus[j])  
  }
  
  # Jointly - Violin plots PER EACH SUPER POPULATION
  for (j in 1:length(l_superpopus)){
    plot_violin_ancestry(df_merged %>% filter(superpopulation %in% l_superpopus[j]),
                         l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, l_superpopus[j])  
  }
  
  # Summary for each locus across all continental groups
  for (j in 1:length(l_superpopus)){
    compute_summary_repeat_per_locus(df_merged %>% filter(superpopulation %in% l_superpopus[j]), 
                                     l_loci[i], 
                                     output_folder)
  }
  
}




