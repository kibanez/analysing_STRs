# Objective: is to see the distribution of repeat-expansions on HTT across 100K
# 1) Definining 3 different threshold cutoffs
# 2) Definining 3 different sub-cohorts or sub-datasets within 85K
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/EHv3.1.2/")

# Load data
dedup_data = read.csv("table_STR_repeat_size_each_row_allele_EHv3.1.2_HTT_CAG_simplified_dedup_270120.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(dedup_data)
# 130126  19




