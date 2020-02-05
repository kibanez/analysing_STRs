# Objective: remove duplicate records from cases-controls files: WHY? There are participants that have been sequenced several times...
# Strategy: we select the latest platekey

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.2 (2019-12-12)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/EHv2.5.5/")

# Load data
merged_data = read.csv("table_STR_repeat_size_each_row_allele_EHv2.5.5_AR_CAG_simplified.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(merged_data)
# 135908  19

# What it happens now, is that `specific disease` are in `specific_disease`, `disease_group`, and `disease_subgroup`
list_spec_disease = merged_data %>% 
  group_by(participant_id) %>% 
  summarise(spec_disease_list = toString(specific_disease)) %>% ungroup() %>% as.data.frame()
dim(list_spec_disease)
# 75095  2

list_disease_group = merged_data %>% 
  group_by(participant_id) %>% 
  summarise(disease_group_list = toString(disease_group)) %>% ungroup() %>% as.data.frame()
dim(list_disease_group)
# 75095  2

list_disease_subgroup = merged_data %>% 
  group_by(participant_id) %>% 
  summarise(disease_subgroup_list = toString(disease_sub_group)) %>% ungroup() %>% as.data.frame()
dim(list_disease_subgroup)
# 75095  2

merged_data = left_join(merged_data,
                        list_spec_disease,
                        by = "participant_id")
dim(merged_data)
# 135908  20

merged_data = left_join(merged_data,
                        list_disease_group,
                        by = "participant_id")
dim(merged_data)
# 135908     21

merged_data = left_join(merged_data,
                        list_disease_subgroup,
                        by = "participant_id")
dim(merged_data)
# 135908     22

# let's remove the other columns now
merged_data = merged_data[,-c(7:9)]
dim(merged_data)
# 135908  19

# There are some participants for which there are several genomes/platekeys 
# we will take/select the latest one
merged_data %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 75095
merged_data %>% select(platekey) %>% unique() %>% pull() %>% length()
# 84917

# First, select duplicated genomes
duplicated_genomes = merged_data %>% 
  group_by(participant_id) %>% 
  filter(n()>2)

duplicated_genomes %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 6483 genomes that have >1 genome
l_duplicated_genomes = unique(duplicated_genomes$participant_id)
length(l_duplicated_genomes)
# 6483

df_pid_platekey = merged_data %>% 
  filter(participant_id %in% l_duplicated_genomes) %>% 
  select(participant_id, platekey)

# From here, select the latest genome
df_pid_platekey = df_pid_platekey %>% 
  group_by(participant_id) %>%
  mutate(latest_platekey = max(platekey)) %>%
  ungroup() %>% 
  as.data.frame()

# select latest genomes
l_latest_dedup_platekeys = unique(df_pid_platekey$latest_platekey)
length(l_latest_dedup_platekeys)
# 6483 (== dedup participant ids)
