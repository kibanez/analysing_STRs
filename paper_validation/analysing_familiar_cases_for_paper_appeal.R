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

# Load main clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V12_and_Pilot_programmes.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(clin_data)
# 2472865  26

# Retrieve the family IDs for all these ~11k PIDs
df_families = clin_data %>%
  filter(participant_id %in% l_neuro) %>%
  select(rare_diseases_family_id, participant_id, affection_status, biological_relationship_to_proband) %>%
  unique()
dim(df_families)
# 11594  4

df_families = df_families %>%
  group_by(participant_id) %>%
  mutate(is_affected = ifelse(affection_status == "Unaffected", FALSE, TRUE)) %>%
  ungroup() %>%
  as.data.frame()

df_families = df_families %>%
  group_by(rare_diseases_family_id, is_affected) %>%
  mutate(num_affected_members = n()) %>%
  ungroup() %>%
  as.data.frame()

