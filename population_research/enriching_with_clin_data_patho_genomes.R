# Script to enrich with clinical data genomes that have an expansion beyond the full-mutation cut-off
# and see whether they can be considered in our cohort
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/feb2021/enriched_with_clin_data/")

# Load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V12_and_Pilot_programmes.tsv", stringsAsFactors = F, header = T, sep= "\t")
dim(clin_data)
# 