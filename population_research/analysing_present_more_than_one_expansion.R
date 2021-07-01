# Objective: compute and analyse how many genomes present an expansion beyond premut and patho thrsehold
# After visual inspection
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# set working dire
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/feb2021/")


