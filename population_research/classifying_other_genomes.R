# Objective: analyse repeat-sizes across other ancestries in the 1Kg genome project
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/1kg/")

# Load data
all_data = read.csv("1000G_2504_high_coverage.sequence.index.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(all_data)
# 2504  22