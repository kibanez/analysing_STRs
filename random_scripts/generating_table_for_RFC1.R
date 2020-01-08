# Objective: generate table with <PLATEKEY> | <VCF_FILE> | <PATH> | <BUILD>
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/CANVAS_RFC1/")

