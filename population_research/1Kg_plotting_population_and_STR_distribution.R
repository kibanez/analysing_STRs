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


# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/1Kg/")

# Load 1Kg population index data
popu_info = read.csv("./integrated_call_samples_v2.20130502.ALL.ped",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)
dim(popu_info)
# 3691  17

# Functions
source("/Users/kibanez/git/analysing_STRs/functions/plot_violin_ancestry.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene_joint_ancestries.R")

# Load EHv3.2.2 STR merged data for each sub-population
df_merged = data.frame()
l_popus = unique(popu_info$Population)
for (i in 1:length(l_popus)){
  popu_aux = paste("~/Documents/STRs/ANALYSIS/population_research/1Kg/data/", l_popus[i] ,sep = "")
  file_aux = paste(paste("merged_", l_popus[i], sep = ""), "_1Kg_samples.tsv", sep = "")
  
  
  df_aux = read.csv(paste(popu_aux, file_aux, sep = "/"),
                    sep  = "\t",
                    stringsAsFactors = F,
                    header = T)
  
  df_aux = df_aux %>% 
    mutate(population = l_popus[i])
  
  df_merged = rbind(df_merged,
                    df_aux)
  
}


df_acb = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/ACB/merged_ACB_1Kg_samples.tsv',
              sep = '\t',
              stringsAsFactors = F,
              header = T)
dim(df_acb)
# 

df_asw = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T)
dim(df_asw)
#

df_beb = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T)
dim(df_beb)

df_ceu = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_ceu)
#

df_clm = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_clm)
#


df_esn = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_esn)
#

df_fin = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_fin)
#

df_gbr = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_gbr)
#

df_gih = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_gih)
#

df_gwd = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_gwd)
#

df_ibs = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_ibs)
#

df_itu = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_itu)
#

df_lwk = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_lwk)
#

df_msl = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_msl)

df_mxl = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_mxl)

df_pel = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_pel)

df_pjl = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_pjl)
#

df_pur = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_pur)

df_stu = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_stu)
#

df_tsi = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_tsi)
#

df_yri = read.csv('~/Documents/STRs/ANALYSIS/population_research/1Kg/data/',
               sep = '\t',
               stringsAsFactors = F,
               header = T)
dim(df_yri)
#


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
# 7290  13

# Population enriched genomes are only GRCh38, we will ignore then GRCh37
output_folder = "./figures/"

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
  
  # Jointly - Violin plots
  plot_violin_ancestry(df_all, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
}