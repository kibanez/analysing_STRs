# Objective: visualise expanded genomes from batch march 2020, across 13 loci, EHv322
# https://docs.google.com/spreadsheets/d/1cuh2rsDkQP3YEHjX6ogLWlgxO3sd0Jfs4zCVdizBQEc/edit#gid=621033475
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.1"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/")

# Load data - all merged
merged_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - merged_all.tsv",
                        stringsAsFactors = F, 
                        header = T,
                        sep = "\t")
dim(merged_table)
# 514  4

# Load MAIN and PILOT ancestry info, to retrieve PC1 and PC2 values, for each genome to be plotted
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      stringsAsFactors = F, 
                      sep = ",",
                      header = T)
dim(popu_table)
# 59464  36

pilot_popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            sep = ",",
                            header = T)
dim(pilot_popu_table)
# 4821  44 

# recode PC names in pilot_popu_table
colnames(pilot_popu_table) = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", colnames(pilot_popu_table)[c(8:44)])

# We cannot Merge Main and Pilot ancestry tables, since PCs have been computed in a diff way

# enrich merged_table with PC1 and PC2 values
list_exp_genomes = unique(merged_table$ID)
length(list_exp_genomes)
# 512

merged_table_main = left_join(merged_table,
                              popu_table %>% filter(ID %in% list_exp_genomes) %>% select(ID, PC1, PC2),
                              by = "ID")
dim(merged_table_main)
# 514  6

merged_table_pilot = left_join(merged_table,
                              pilot_popu_table %>% filter(ID %in% list_exp_genomes) %>% select(ID, PC1, PC2),
                              by = "ID")
dim(merged_table_pilot)
# 514  6

png("figures/expanded_genomes_MAIN.png")
ggplot(data=merged_table_main %>% filter(!is.na(merged.superpopu)), 
       aes(x=PC2, y=PC1, colour = merged.superpopu)) +
  geom_point() +
  xlab("PC2") +
  ylab("PC1") +
  guides(fill = FALSE)
dev.off()

png("figures/expanded_genomes_PILOT.png")
ggplot(data=merged_table_pilot %>% filter(!is.na(merged.superpopu), !is.na(PC1)), 
       aes(x=PC2, y=PC1, colour = merged.superpopu)) +
  geom_point() +
  xlab("PC2") +
  ylab("PC1") +
  guides(fill = FALSE)
dev.off()

