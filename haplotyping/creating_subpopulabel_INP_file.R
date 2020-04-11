# Objective: create subpopuation labels INP file for fastPHASE input
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/AR/")
