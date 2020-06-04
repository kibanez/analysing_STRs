# Objective: visualise expanded genomes from batch march 2020, across 13 loci, EHv322
# https://docs.google.com/spreadsheets/d/1cuh2rsDkQP3YEHjX6ogLWlgxO3sd0Jfs4zCVdizBQEc/edit#gid=621033475
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.1"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/")

# Load data - all merged
merged_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - merged_all.tsv",
                        stringsAsFactors = F, 
                        header = T,
                        sep = "\t")
dim(merged_table)
# 514  4