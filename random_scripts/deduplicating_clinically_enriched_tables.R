# Objective: remove duplicate records from cases-controls files: WHY? There are participants that have been sequenced several times...
# Strategy: we select the latest platekey
# In particular, this script has been implemented because of HTT and Henrietta's project
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/EHv2.5.5/")

# Load data
merged_data = read.csv("table_STR_repeat_size_each_row_allele_EHv2.5.5_HTT_CAG_simplified.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(merged_data)
# 176896     19


# There are some participants for which there are several genomes/platekeys 
# we will take/select the latest one

merged_data %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 76148
merged_data %>% select(plate_key.x) %>% unique() %>% pull() %>% length()
# 86243

# First, select duplicated genomes
duplicated_genomes = merged_data %>% 
  group_by(participant_id) %>% 
  filter(n()>3)

duplicated_genomes %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 11283 genomes that have >1 genome
l_duplicated_genomes = unique(duplicated_genomes$participant_id)
length(l_duplicated_genomes)
# 11283

df_pid_platekey = merged_data %>% filter(participant_id %in% l_duplicated_genomes) %>% select(participant_id, plate_key.x)

# From here, select the latest genome
df_pid_platekey = df_pid_platekey %>% 
  group_by(participant_id) %>%
  mutate(latest_platekey = max(plate_key.x)) %>%
  ungroup() %>% 
  as.data.frame()

# select latest genomes
l_latest_dedup_platekeys = unique(df_pid_platekey$latest_platekey)
length(l_latest_dedup_platekeys)
# 11283 (== dedup participant ids)

# PART 1 - select pid not duplicated (to merge with the rest afterwards)
merged_data_dedup = merged_data %>% filter(!participant_id %in% l_duplicated_genomes)
dim(merged_data_dedup)
# 129730  19

# PART 2 - include the genomes recovered from the duplicated genomes
merged_data_dedup = rbind(merged_data_dedup,
                          merged_data %>% filter(plate_key.x %in% l_latest_dedup_platekeys))
dim(merged_data_dedup)
# 156590  19


