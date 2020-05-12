# Objective: define algorithmia behind annotated case-control TSV file when selecting repeat-motifs that are enriched in cases rather than in controls
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.9.0/case-control/analysis/")

# Loading data
cc_data = read.csv("./test_ALS/output/casecontrol_test_ALS_locus-based_annotated.tsv",
                   stringsAsFactors = F,
                   header = T,
                   sep = "\t")
dim(cc_data)
# 8832  9

# List of cases
l_cases = read.table("./test_ALS/input/main_6_cases.txt",
                     stringsAsFactors = F)
l_cases = l_cases$V1
length(l_cases)
# 6

# List of controls
l_controls = read.table("./test_ALS/input/main_11355_controls.txt", stringsAsFactors = F)
l_controls = l_controls$V1
length(l_controls)
# 11355

# Filtering criteria: pvalue <= 0.05 
cc_data_filtered = cc_data %>%
  filter(pvalue <= 0.05)
dim(cc_data_filtered)
# 186  9

# Filtering those genes or repeat-motifs that are significantly enriched in cases
# We need to analyse `counts` field
# Create cases and controls columns, and put there the % of how many are there??


