# Objective: retrieve the genomes (platekeys) for which EH (both v2.5.5 and v3.1.2) have estimated an insertion-STR in RFC1 or Canvas
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/CANVAS_RFC1/")
