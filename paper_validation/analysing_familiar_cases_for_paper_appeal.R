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

# Load merged PID-FID clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/merged_RE_releases_and_Pilot_RD_and_Cancer_PID_FID_platekey_up_to_RE_V12.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(clin_data)
# 175781  6

l_families_subcohort = clin_data %>% filter(participant_id %in% l_neuro) %>% select(rare_diseases_family_id) %>% unique() %>% pull()
length(l_families_subcohort)
# 10417

# Load pedigree table (MAIN)
pedigree_table = read.csv("~/Documents/STRs/clinical_data/clinical_data/raw/rare_diseases_pedigree_2021-07-15_12-35-12.tsv",
                          stringsAsFactors = F,
                          header = T,
                          sep = "\t")
dim(pedigree_table)
# 34459  4

# Load pedigree-member table (MAIN)
pedigree_member_table = read.csv("~/Documents/STRs/clinical_data/clinical_data/raw/rare_diseases_pedigree_member_2021-05-24_13-30-57.tsv",
                                 stringsAsFactors = F,
                                 header = T,
                                 sep = "\t")
dim(pedigree_member_table)
# 217704  35 

pedigree_merged = left_join(pedigree_table,
                            pedigree_member_table %>% select(participant_id, rare_diseases_pedigree_sk, rare_diseases_family_sk, proband, life_status, affection_status),
                            by = "rare_diseases_family_sk")
dim(pedigree_merged)
# 217779  9



# Retrieve the family IDs for all these ~11k PIDs
df_families = clin_data %>%
  filter(rare_diseases_family_id %in% l_families_subcohort) %>%
  select(rare_diseases_family_id, participant_id, affection_status, biological_relationship_to_proband) %>%
  unique()
dim(df_families)
# 24749  4

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

# Compute how many families we've got with more than 1 affected members
df_families %>% filter(num_affected_members > 1) %>% select(rare_diseases_family_id) %>% unique() %>% pull() %>% length()
# 6145cl



