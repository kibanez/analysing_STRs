# Objective: there is a list of genomes that have a pathogenic expansion on FXN, but only in one allele. 
# we want to analyse whether they also have a potentially pathogenic variant in FXN
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/FXN_monoallelic/")

