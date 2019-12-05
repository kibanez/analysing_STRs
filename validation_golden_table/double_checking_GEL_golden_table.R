# Objective: we have filled in GEL genomes in golden validation table. We now want to check whether classification_a1 and classification_a2 for EHv2 and EHv3 are correct
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# set working directory
setwd("~/Documents/STRs/VALIDATION/")

# load latest validation golden table
val_data = read.csv("EHv255_EHv312_validation_cohort_GEL.tsv",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)
dim(val_data)
# 638  28

