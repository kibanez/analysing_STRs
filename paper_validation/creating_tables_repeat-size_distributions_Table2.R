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
dim(atn1_table)
# 23864  19

length(unique(atn1_table$platekey))
# 12994
length(unique(atn1_table$participant_id))
# 12994

write.table(atn1_table,
            "./case-controls_ATN1_12994_pids_removingUnique.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# ATXN2
atxn2_table = read.csv("~/Documents/STRs/ANALYSIS/cases_controls/batch_march/EHv255/table_STR_repeat_size_each_row_allele_EHv2.5.5_ATXN2_CAG_simplified.tsv",
                       stringsAsFactors = F,
                       sep = "\t",
                       header = T)
dim(atxn2_table)
# 188714  19

atxn2_table = atxn2_table %>%
  filter(participant_id %in% l_pid)
dim(atxn2_table)
# 16188 19

length(unique(atxn2_table$platekey))
# 12993
length(unique(atxn2_table$participant_id))
# 12993

write.table(atxn2_table,
            "./case-controls_ATXN2_12993_pids_removingUnique.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)


# ATXN7
atxn7_table = read.csv("~/Documents/STRs/ANALYSIS/cases_controls/batch_march/EHv255/table_STR_repeat_size_each_row_allele_EHv2.5.5_ATXN7_CAG_simplified.tsv",
                       stringsAsFactors = F,
                       sep = "\t",
                       header = T)
dim(atxn7_table)
# 188680  19

atxn7_table = atxn7_table %>%
  filter(participant_id %in% l_pid)
dim(atxn7_table)
# 17987  19

length(unique(atxn7_table$participant_id))
# 12993
length(unique(atxn7_table$platekey))
# 12993

write.table(atxn7_table,
            "./case-controls_ATXN7_12993_pids_removingUnique.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)



# HTT
htt_table = read.csv("~/Documents/STRs/ANALYSIS/cases_controls/batch_march/EHv255/table_STR_repeat_size_each_row_allele_EHv2.5.5_HTT_CAG_simplified.tsv",
                     stringsAsFactors = F,
                     sep = "\t",
                     header = T)
dim(htt_table)
# 188718  19

htt_table = htt_table %>%
  filter(participant_id %in% l_pid)
dim(htt_table)
# 24119 19

length(unique(htt_table$participant_id))
# 12994
length(unique(htt_table$platekey))
# 12994

write.table(htt_table,
            "./case-controls_HTT_12994_pids_removingUnique.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

