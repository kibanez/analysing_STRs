# Objective: from the work we have done by inspecting visually all pileups, compute the carrier ratio for each locus
# unrelated, probands unrelated, probands not neuro
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/")
