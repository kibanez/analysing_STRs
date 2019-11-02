# Objective: Complete our validation table together with the max CI value estimated by EH for each allele

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr)

# set the working directory
setwd("~/Documents/STRs/clinical_validation_paper/plots/")
