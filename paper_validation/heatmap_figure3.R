# Objective: Create a heatmap with the number of confirmed repeat-expansions across loci and diseases
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"

# set the working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/LANCET/APPEAL/")

# Load raw data with the confirmed RE across diseases and genes
table3 = read.csv("./template_for_figure3.tsv",
                  stringsAsFactors = F, 
                  header = T,
                  sep = "\t")
dim(table3)
# 15  14