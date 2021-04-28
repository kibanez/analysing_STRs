# Objective: compute repeat-size concordance between EH and PCR, just keeping 1 gene for those duplicated or repeated genomes
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"

# set the working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

concordance = read.csv("./concordance_repeats.tsv", stringsAsFactors = F, sep = "\t", header = T)
dim(concordance)
# 902  11

platekey_duplicated = unique(concordance$platekey[which(duplicated(concordance$platekey))])
length(platekey_duplicated)
# 195

# There are a total of 195 repeat platekeys - for which we've got PCR lengths for more than one gene