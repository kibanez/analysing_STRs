# Objective; calculate the total number of unique participants across A,B,C, and D subtables
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"

# set environemnt
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/subtables/")

# Load final complete subtables
table_a_main = read.csv("./TableA_main_including_ultrarare_same_cutoff_for_tableA.tsv",
                        stringsAsFactors = F,
                        sep = "\t",
                        header = T)
table_a_pilot = read.csv("./TableA_pilot_new_cutoffs.tsv",
                         stringsAsFactors = F,
                         header = T,
                         sep = "\t")

table_b_main = read.csv("./TableB_main.tsv",
                        stringsAsFactors = F,
                        header = T,
                        sep = "\t")

table_c_main = read.csv("./TableC_main.tsv",
                        stringsAsFactors = F,
                        header = T,
                        sep = "\t")

