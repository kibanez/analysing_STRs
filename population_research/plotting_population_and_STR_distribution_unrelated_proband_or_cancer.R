# Objective: take advantage from the work GEL bioinfo research group has done estimating ancestry for ~59,356 genomes
# and analyse and study the STR distribution across them
# Only in a subset with unrelated, affected/proband or cancer genomes
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
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/")

# Load population data
popu_table_enriched = read.csv("./population_info_enriched_59356_by_031019.tsv",
                               header = T,
                               sep = "\t",
                               stringsAsFactors = F)
dim(popu_table_enriched)
# 59356  20

# Load thresholds
# STR annotation, threshold including the largest normal and the smallest pathogenic sizes reported
gene_annotation_normal = '/Users/kibanez/git/analysing_STRs/threshold_largest_normal_reported_research.txt'
gene_data_normal = read.table(gene_annotation_normal, stringsAsFactors=F, header = T)

gene_annotation_pathogenic = '/Users/kibanez/git/analysing_STRs/threshold_smallest_pathogenic_reported_research.txt'
gene_data_pathogenic = read.table(gene_annotation_pathogenic, stringsAsFactors=F, header = T)


# Functions
source("/Users/kibanez/git/analysing_STRs/functions/plot_violin_ancestry.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene_joint_ancestries.R")
source("/Users/kibanez/git/analysing_STRs/functions/compute_summary_repeat_per_locus.R")

# Load EH STR output data - merged TSV files
# Research ~80K genomes, EH-v3.1.2- November 2019 -- but the ones that are ALSO included in the population info table
df_asi = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/ASI/merged/merged_population_genomes_unrelated_probands_and_cancer_2678_avg_EHv3.1.2_ASI.tsv',
              sep = '\t',
              stringsAsFactors = F,
              header = T)
dim(df_asi)
# 1430  12

df_eas = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/EAS/merged/merged_population_genomes_unrelated_probands_and_cancer_227_avg_EHv3.1.2_EAS',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_eas)
# 723  12

df_afr = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/AFR/merged/merged_population_genomes_unrelated_probands_and_cancer_1136_avg_EHv3.1.2_AFR.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_afr)
# 1170  12

df_amr = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/AMR/merged/merged_population_genomes_unrelated_probands_and_cancer_384_avg_EHv3.1.2_AMR.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_amr)
# 965  12

df_eur = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/EUR/merged/merged_population_genomes_unrelated_probands_and_cancer_26033_avg_EHv3.1.2_EUR.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_eur)
# 2282  12


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
# 6570  13

# Population enriched genomes are only GRCh38, we will ignore then GRCh37
output_folder = "./figures_unrelated_affected_or_cancer/"

l_loci = sort(unique(df_all$gene))
for (i in 1:length(l_loci)){
  # Each locus - Individually
  plot_gene(df_afr, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "AFR")
  plot_gene(df_amr, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "AMR")
  plot_gene(df_eur, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "EUR")
  plot_gene(df_eas, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "EAS")
  plot_gene(df_asi, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "ASI")
  
  # Jointly - distribution
  plot_gene_joint_ancestries(df_all, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
  
  # Jointly - Violing plots
  plot_violin_ancestry(df_all, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
  
  # Summary for each locus across all continental groups
  compute_summary_repeat_per_locus(df_all, l_loci[i], output_folder)
  
}
