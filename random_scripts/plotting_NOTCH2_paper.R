# Objective: plot percentage (not frequency) of each repeat-size in 
# - individuals recruited under `neurology`
# - rest of individuals
# - NIID
# - inclusions
# last both sent by Zhongbo
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.2 (2019-12-12)"

# libraries
library(dplyr)
library(ggplot2)

# Set working directory
working_dir="~/Documents/STRs/PAPERS/NIID/"
setwd(working_dir)

# Load data
df_zhongbo = read.csv("data/repeat_sizes_NIID_and_inclusions_for_Kristina.csv",
                      sep = ",",
                      header = T,
                      stringsAsFactors = F)
dim(df_zhongbo)
# 30  6

# Load data corresponding to 56K genomes (to align with other figures) we do have for populations
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/population_and_super-population_definitions_across_59352_WGS_REv9_271119.tsv",
                               header = T,
                               sep = "\t",
                               stringsAsFactors = F)
dim(popu_table)
# 59356  21

l_genomes = unique(popu_table$platekey)
length(l_genomes)
# 59356



joint_plot = ggplot(df_gene_barplot, aes(x = number_repeats, y = af, group = population, color = population)) + 
  geom_line() + 
  ylab("Allele frequency") + 
  xlab("Repeat sizes (repeat units)") + 
  ggtitle(gene_name) + 
  geom_vline(xintercept = threshold_normal, colour = 'blue', lty = 2) + 
  geom_vline(xintercept = threshold_pathogenic, colour = 'red', lty = 2) + 
  coord_cartesian(xlim = c(min_value,max_value))
