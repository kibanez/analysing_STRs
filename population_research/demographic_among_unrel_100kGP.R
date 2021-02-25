# Objective: study and retrieve numbers re demographic within 100kGP
# unrelated genomes, filtering out 125bp genomes
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(tidyverse)

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/")

# Load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V10_and_Pilot_programmes.tsv",
                     stringsAsFactors = F, 
                     header = T,
                     sep = "\t")
dim(clin_data)
# 2101385  24

# Load unrel genomes
l_unrel = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/l_unrelated_55603_genomes_batch2.txt",
                     stringsAsFactors = F)
l_unrel = l_unrel$V1
length(l_unrel)
# 55603
