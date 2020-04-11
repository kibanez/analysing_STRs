# Objective: select the best K, number of clusters, when running fastPHASE
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/AR/")

