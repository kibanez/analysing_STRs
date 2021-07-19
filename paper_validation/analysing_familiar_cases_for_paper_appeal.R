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
                            pedigree_member_table %>% select(participant_id, rare_diseases_family_sk, rare_diseases_pedigree_member_sk, proband, life_status, affection_status),
                            by = "rare_diseases_family_sk")
dim(pedigree_merged)
# 217779  9

pedigree_merged = pedigree_merged %>%
  group_by(rare_diseases_pedigree_sk) %>%
  mutate(n_affected_members = sum(affection_status == "Affected")) %>%
  ungroup() %>%
  as.data.frame()

# Retrieve the family IDs for all these ~11k PIDs
df_families = pedigree_merged %>%
  filter(rare_diseases_family_id %in% l_families_subcohort) %>%
  unique()
dim(df_families)
# 67006  10

# Define whether a family has family-history (num_affected_members > 1)
df_families = df_families %>%
  group_by(rare_diseases_family_id) %>%
  mutate(is_fam = ifelse(n_affected_members > 1, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame() %>%
  unique()

# Compute how many families we've got with more than 1 affected members
df_families %>% filter(n_affected_members > 1) %>% select(rare_diseases_family_id) %>% unique() %>% pull() %>% length()
# 1847 (Main programme: 17.73%)
df_families %>% filter(is_fam) %>% select(rare_diseases_family_id) %>% unique() %>% pull() %>% length()
# 1847

l_pids_familials = df_families %>% filter(num_affected_members > 1) %>% select(rare_diseases_family_id.x) %>% unique() %>% pull()
write.table(l_pids_familials, "~/Documents/STRs/PAPERS/VALIDATION_PAPER/LANCET/APPEAL/list_5174_families_at_least_2_affected.tsv",
            quote = F, col.names = F, row.names = F)

# Load pilot clinical data
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           stringsAsFactors = F,
                           header = T,
                           sep = "\t")
dim(pilot_clin_data)
# 4974  10

l_pilot_fam = intersect(l_families_subcohort,
                        pilot_clin_data$gelFamilyId.x)
length(l_pilot_fam)
# 557

# Enrich pilot clin data with number of affected members
pilot_clin_data_subcohort = pilot_clin_data %>%
  filter(gelFamilyId.x %in% l_pilot_fam) %>%
  group_by(gelFamilyId.x) %>%
  mutate(num_affected_members = sum(disease_status %in% c("Affected", "AffectedSame"))) %>%
  ungroup() %>%
  as.data.frame() %>%
  unique()
dim(pilot_clin_data_subcohort)
# 1308

# How many pilot families have at least 1 other member in the family affected?
pilot_clin_data_subcohort %>% filter(num_affected_members > 1) %>% select(gelFamilyId.x) %>% unique() %>% pull() %>% length()
# 125

write.table(df_families, "table_family_history_for_11k_PIDs_MAIN_PROGRAMME.tsv", 
            quote = F, row.names = F, col.names = T, sep= "\t")
write.table(pilot_clin_data_subcohort, "table_family_history_for_11k_PIDs_PILOT_PROGRAMME.tsv", 
            quote = F, row.names = F, col.names = T, sep= "\t")

l_pids_familials_pilot = pilot_clin_data_subcohort %>% filter(num_affected_members > 1) %>% select(gelFamilyId.x) %>% unique() %>% pull() 
write.table(l_pids_familials_pilot,
            "~/Documents/STRs/PAPERS/VALIDATION_PAPER/LANCET/APPEAL/list_125_families_PILOT_at_least_2_affected.tsv",
            quote = F, col.names = F, row.names = F)

# Total number: MAIN + PILOT
#1847 + 125 =  1972 (19%)

# What about all the cohort?
pedigree_merged = pedigree_merged %>%
  group_by(rare_diseases_pedigree_sk) %>%
  mutate(num_affected_members = sum(affection_status == "Affected")) %>%
  ungroup() %>%
  as.data.frame() %>%
  unique()

pedigree_merged %>% filter(num_affected_members > 1) %>% select(rare_diseases_family_id) %>% unique() %>% pull() %>% length()
# 8397

length(unique(pedigree_merged$rare_diseases_family_id))
# 34453 (24.4%)
