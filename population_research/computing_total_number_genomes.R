# Script to compute total number of genomes in the 100kGP
# We filter out genomes sequenced in read-length 125bp
# All vs not considering neuro
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/")

l_125 = read.table("./list_genomes_125bp_100kGP.tsv", stringsAsFactors = F)
l_125 = l_125$V1

l_unrel = read.table("./l_unrelated_55603_genomes_batch2.txt", stringsAsFactors = F)
l_unrel = l_unrel$V1

l_unrel_not125 = unique(setdiff(l_unrel, l_125))
length(l_unrel_not125)
# 54437