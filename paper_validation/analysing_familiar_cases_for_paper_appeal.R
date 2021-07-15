# Objective: analyse how many family-cases we've got inside ~11k participants
# a familiar case means when anyone within the family-history (regardless of whether it's been recruited or not) is affected
# for appeal work - july 2021
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

# load dataset
l_neuro = read.table("./list_11631_PIDs.txt", stringsAsFactors = F)
l_neuro = l_neuro$V1
length(l_neuro)
# 11631

# For MAIN programme we need to use `rare_disease_pedigree_member` table - `affection_status` column/parameter
# For PILOT programme we need to use `pedigree` table - `status` column/parameter
main_clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/raw/rare_diseases_pedigree_member_2021-05-24_13-30-57.tsv",
                          stringsAsFactors = F,
                          header = T,
                          sep = "\t")
dim(main_clin_data)
# 217704 35

pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/data_freeze_Pilot_LK_RESEARCH/pedigree.csv",
                           stringsAsFactors = F,
                           header = T)
dim(pilot_clin_data)
# 17258  34

# Create df with PID and the affection status in main or pilot
df_affected_main = main_clin_data %>% 
  filter(participant_id %in% l_neuro) %>%
  select(participant_id, affection_status) %>%
  unique()
dim(df_affected_main)
# 10615  2

df_affected_pilot = pilot_clin_data %>%
  filter(gelId %in% l_neuro) %>%
  select(gelId, affectionStatus) %>%
  unique()
dim(df_affected_pilot)
# 650  2
colnames(df_affected_pilot) = colnames(df_affected_main)

df_affected_merged = unique(rbind(df_affected_main,
                                  df_affected_pilot))
dim(df_affected_merged)
# 11265   2