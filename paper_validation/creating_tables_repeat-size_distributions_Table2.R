# Objective: count the total number of participants, %gender. %ethnicity 
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"

# set the working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

# Load merged data
l_pid = read.table("./list_13331_PIDs_table2.txt", stringsAsFactors = F)
l_pid = l_pid$V1
length(l_pid)
# 13331

# Strategy: take EHv255 case-control tables from latest batch march 2020
# ATN1
atn1_table = read.csv("~/Documents/STRs/ANALYSIS/cases_controls/batch_march/EHv255/table_STR_repeat_size_each_row_allele_EHv2.5.5_ATN1_CAG_simplified.tsv",
                    stringsAsFactors = F,
                    sep = "\t",
                    header = T)
dim(atn1_table)
# 188718  19

atn1_table = atn1_table %>%
  filter(participant_id %in% l_pid)
atn1_table  = unique(atn1_table)
dim(atn1_table)
# 23864  19

length(unique(atn1_table$platekey))
# 12994
length(unique(atn1_table$participant_id))
# 12994

write.table(atn1_table,
            "./case-controls_ATN1_12994_pids.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
