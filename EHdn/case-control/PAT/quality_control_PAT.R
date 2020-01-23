# Objective: to see age distribution across case and control datasets
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/case-control/PAT/")

# pilot clinical data
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           sep = "\t",
                           stringsAsFactors = F,
                           header = T)
dim(pilot_clin_data)
# 4974  10

# main clinical data
main_clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T)
dim(main_clin_data)
# 1124633  28
