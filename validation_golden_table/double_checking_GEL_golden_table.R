# Objective: we have filled in GEL genomes in golden validation table. We now want to check whether classification_a1 and classification_a2 for EHv2 and EHv3 are correct
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
