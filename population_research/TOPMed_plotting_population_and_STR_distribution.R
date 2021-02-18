# Objective: analyse the repeat-sizes across sub-population and super-population within TOPMed
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 2020-02-29)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"
library(scatterplot3d); packageDescription("scatterplot3d", fields = "Version") # 0.3-41
library(ggpubr); packageDescription("ggpubr", fields = "Version") # 0.2.3
library(tidyverse)

# Set environment
setwd("~/Documents/STRs/ANALYSIS/population_research/TOPMed/After_QC/AllUnrelatedSamples/")

# Load thresholds
# STR annotation, threshold including the largest normal and the smallest pathogenic sizes reported
gene_annotation_normal = '/Users/kibanez/git/analysing_STRs/threshold_largest_normal_reported_research.txt'
gene_data_normal = read.table(gene_annotation_normal, stringsAsFactors=F, header = T)

gene_annotation_pathogenic = '/Users/kibanez/git/analysing_STRs/threshold_smallest_pathogenic_reported_research.txt'
gene_data_pathogenic = read.table(gene_annotation_pathogenic, stringsAsFactors=F, header = T)

# Functions
source("/Users/kibanez/git/analysing_STRs/functions/plot_violin_ancestry_gnomAD.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene_joint_ancestries_gnomAD.R")
source("/Users/kibanez/git/analysing_STRs/functions/compute_summary_repeat_per_locus.R")

# Load total number of genomes per ancestry
total_num_genomes = read.csv("./number_genomes_per_ancestry.tsv",
                             stringsAsFactors = F, 
                             header = T,
                             sep = "\t")

# Load EHv3.2.2 STR merged data for each sub-population
l_popus = c("AFR", "AMR", "EAS", "EUR", "SAS")
df_merged = data.frame()
for (i in 1:length(l_popus)){
  popu_aux = paste(l_popus[i],"AfterPlotQC_13LociAlleleFreq_WOSampleIds.tsv" ,sep = "_")
  
  num_genomes_ancestry = total_num_genomes %>%
    filter(superpopu %in% l_popus[i]) %>%
    select(list_vcf) %>%
    pull()
  
  df_aux = read.csv(popu_aux,
                    stringsAsFactors = F, 
                    header = T,
                    sep = "\t")
  df_aux$superpopu = rep(l_popus[i], length(df_aux$chr))
  
  df_aux$AF = df_aux$num_samples / (2*num_genomes_ancestry)

  df_merged = rbind(df_merged,
                    df_aux)
}

dim(df_merged)
# 2263  12

output_folder = "./figures/"
l_loci = sort(unique(df_merged$gene))
for (i in 1:length(l_loci)){
  colnames(df_merged)[12] = "superpopu"
  # Specifying sub-population  
  for (j in 1:length(l_popus)){
    # Each locus - Individually
    plot_gene(df_merged %>% filter(superpopu %in% l_popus[j]), l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", l_popus[j])  
  }
  
  # Jointly - distribution
  plot_gene_joint_ancestries_gnomAD(df_merged, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
  
  # Jointly - Violin plots
  plot_violin_ancestry_gnomAD(df_merged, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
  
  # Summary for each locus across all continental groups
  colnames(df_merged)[12] = "population"
  compute_summary_repeat_per_locus(df_merged, l_loci[i], output_folder)
}

