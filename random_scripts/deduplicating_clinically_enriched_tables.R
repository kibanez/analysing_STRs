# Objective: remove duplicate records from cases-controls files
# In particular, this script has been implemented because of HTT and Henrietta's project
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/EHv2.5.5/")

# Load data
merged_data = read.csv("table_STR_repeat_size_each_row_allele_EHv2.5.5_HTT_CAG_simplified.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(merged_data)
# 176896     19


