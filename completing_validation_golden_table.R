# Objective: Complete our validation table together with the max CI value estimated by EH for each allele
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr)

# set the working directory
setwd("~/Documents/STRs/VALIDATION/")

# Load validation golden table data
val_data = read.csv("./STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez.tsv",
                    header = T,
                    sep = "\t",
                    stringsAsFactors = F)
dim(val_data)
# 639 17

# Load Pilot merged table
merged_maxCI_table_pilot = read.csv("./raw_data/pilot_validation/merged/merged_validation_pilot.tsv",
                              sep = "\t",
                              header = T,
                              stringsAsFactors = F)
dim(merged_maxCI_table_pilot)
# 292 11

# Load 
