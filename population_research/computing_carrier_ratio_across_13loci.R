# Objective: from the work we have done by inspecting visually all pileups, compute the carrier ratio for each locus
# unrelated, probands unrelated, probands not neuro
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/")

# Load clinical data for 90K pids across 13 loci
clin_data = read.csv("~/Documents/STRs/ANALYSIS/cases_controls/batch_march/EHv322/table_13_loci_across_90863_genomes_each_row_allele_EHv3.2.2_90863platekeys_88827pids.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(clin_data)
# 91857  46

