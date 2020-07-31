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