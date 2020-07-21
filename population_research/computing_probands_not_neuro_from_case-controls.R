# Objective: count how many probands we do have
# Count how many probands MINUS neuro we do have
# For EHv3.2.2
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/batch_march/EHv322/")
