# Objective: from EHv2.5.5 and EHv3.1.2 tables re STR calls, plot corresponding ROC curves

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# set the working directory
setwd("~/Documents/STRs/VALIDATION/raw_data/")

