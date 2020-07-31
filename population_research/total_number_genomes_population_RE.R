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
