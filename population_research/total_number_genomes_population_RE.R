# Objective: Keep the track of all genomes included in population analyses
# Batch 1 - Oct/Nov 2019
# Batch 2 - August 2020
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
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/")

# Batch1
batch1 = read.csv("GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                  stringsAsFactors = F,
                  header = T)
dim(batch1)
# 59464  36

# Batch2 
batch2 = read.csv("batch2/aggV2_M30K_60KupscaledPCs_R9_08062020.tsv",
                  stringsAsFactors = F,
                  header = T,
                  sep = " ")
dim(batch2)
# 78388  21

# List of unique genomes in both batch 1 and 2
l_unique_genomes = unique(c(unique(batch1$ID),
                     unique(batch2$plate_key)))
length(l_unique_genomes)
# 79849

# How many extra new genomes in batch2?
length(unique(setdiff(unique(batch2$plate_key),
                      unique(batch1$ID))))
# 20,385

# How many genomes in batch1 that are NOT in batch2?
length(unique(setdiff(unique(batch1$ID),
                      unique(batch2$plate_key))))
# 1,461

write.table(l_unique_genomes,
            "list_79849_unique_genomes_batch1_batch2.txt",
            quote = F,
            row.names = F,
            col.names = F)
