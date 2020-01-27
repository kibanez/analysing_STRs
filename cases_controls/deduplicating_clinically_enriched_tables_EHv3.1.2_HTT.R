# Objective: remove duplicate records from cases-controls files: WHY? There are participants that have been sequenced several times...
# Strategy: we select the latest platekey

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


# What it happens now, is that `specific disease` are in `specific_disease`, `disease_group`, and `disease_subgroup`
list_spec_disease = merged_data %>% 
  group_by(participant_id) %>% 
  summarise(spec_disease_list = toString(specific_disease)) %>% ungroup() %>% as.data.frame()
dim(list_spec_disease)
# 76148  2

list_disease_group = merged_data %>% 
  group_by(participant_id) %>% 
  summarise(disease_group_list = toString(disease_group)) %>% ungroup() %>% as.data.frame()
dim(list_disease_group)
# 76148  2

list_disease_subgroup = merged_data %>% 
  group_by(participant_id) %>% 
  summarise(disease_subgroup_list = toString(disease_sub_group)) %>% ungroup() %>% as.data.frame()
dim(list_disease_subgroup)
# 76148  2

merged_data = left_join(merged_data,
                      list_spec_disease,
                      by = "participant_id")
dim(merged_data)
# 176896     20

merged_data = left_join(merged_data,
                        list_disease_group,
                        by = "participant_id")
dim(merged_data)
# 176896     21

merged_data = left_join(merged_data,
                        list_disease_subgroup,
                        by = "participant_id")
dim(merged_data)
# 176896     22

# let's remove the other columns now
merged_data = merged_data[,-c(7:9)]
dim(merged_data)
# 176896  19


# There are some participants for which there are several genomes/platekeys 
# we will take/select the latest one
merged_data %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 76148
merged_data %>% select(plate_key.x) %>% unique() %>% pull() %>% length()
# 86243

# First, select duplicated genomes
duplicated_genomes = merged_data %>% 
  group_by(participant_id) %>% 
  filter(n()>2)

duplicated_genomes %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 11283 genomes that have >1 genome
l_duplicated_genomes = unique(duplicated_genomes$participant_id)
length(l_duplicated_genomes)
# 11283

df_pid_platekey = merged_data %>% 
  filter(participant_id %in% l_duplicated_genomes) %>% 
  select(participant_id, plate_key.x)

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
                          merged_data %>% filter(plate_key.x %in% l_latest_dedup_platekeys, genome_build %in% "GRCh38"))
dim(merged_data_dedup)
# 134684  19

# so far, we have removed duplicates having more than 1 genome
# QUALITY CHECK - check whether there are pids with more than 3 rows
merged_data_dedup %>% 
  group_by(participant_id) %>% 
  filter(n()>2) %>%
  dim()
# 1886  19

l_duplicated_genomes2 = merged_data_dedup %>% 
  group_by(participant_id) %>% 
  filter(n()>3) %>%
  select(participant_id) %>%
  pull() %>%
  unique()
length(l_duplicated_genomes2)
# 458

# Take the largets repeat-size
df_pid_platekey_repeatsize = merged_data_dedup %>% 
  filter(participant_id %in% l_duplicated_genomes2) %>% 
  select(participant_id, plate_key.x, repeat_size)

df_pid_platekey_repeatsize = df_pid_platekey_repeatsize %>% 
  group_by(participant_id) %>%
  mutate(largest_repeat = max(repeat_size)) %>%
  ungroup() %>% 
  as.data.frame()

# remove from merged_data_dedup those genomes with smaller repeat-size

# And now, only include here the ones with largest repeat-size
df_aux_recover = merged_data_dedup %>% filter(participant_id %in% l_duplicated_genomes2)

# First, remove all genomes in l_genomes_duplciated2
merged_data_dedup = merged_data_dedup %>% filter(!participant_id %in% l_duplicated_genomes2)
dim(merged_data_dedup)
# 132798  19


df_aux_recover = df_aux_recover %>% group_by(participant_id) %>% mutate(large_repeat = max(repeat_size)) %>% ungroup() %>% as.data.frame()
index_to_keep = which(df_aux_recover$large_repeat == df_aux_recover$repeat_size)
df_aux_recover = df_aux_recover[index_to_keep,]
dim(df_aux_recover)
# 1097  20

# remove last column
df_aux_recover = df_aux_recover[,-20]
dim(df_aux_recover)
# 1097  19

# There are some genomes in GRCh37 and GRCh38, let's take GRCh38 info
l_genomes_both_builds = df_aux_recover %>%
  group_by(participant_id) %>%
  filter(n()>2) %>%
  ungroup() %>%
  select(participant_id) %>%
  unique() %>%
  pull()

df_aux_recover = rbind(df_aux_recover %>% 
                         filter(!participant_id %in% l_genomes_both_builds),
                       df_aux_recover %>%
                         filter(participant_id %in% l_genomes_both_builds, genome_build %in% "GRCh38"))

dim(df_aux_recover)
# 1097  19

merged_data_dedup = rbind(merged_data_dedup,
                          df_aux_recover)
dim(merged_data_dedup)
#  133895  19

# QUALITY CHECK
merged_data_dedup %>% 
  group_by(participant_id) %>% 
  filter(n()>2) %>%
  dim()
# 373  19

aver = merged_data_dedup %>% 
  group_by(participant_id) %>% 
  filter(n()>2) %>%
  ungroup() %>%
  as.data.frame()
length(unique(aver$participant_id))
# 96
# there are 96 participant-genome rows with too many rows...let's simplify them

l_parti_too_many_rows = unique(aver$participant_id)
merged_data_dedup_final = merged_data_dedup %>% 
  filter(!participant_id %in% l_parti_too_many_rows)

aver = unique(aver)
dim(aver)
# 96  19
# repeat x2 each row
aver = aver[rep(seq_len(nrow(aver)), each = 2), ]
dim(aver)
# 192  19

merged_data_dedup_final = rbind(merged_data_dedup_final,
                                aver)

dim(merged_data_dedup_final)
# 133714  19

# How many genomes?
length(unique(merged_data_dedup_final$plate_key.x))
# 66857

# how many participant ids?
length(unique(merged_data_dedup_final$participant_id))
# 66867

# how many alleles? not unique
length(merged_data_dedup_final$repeat_size)
# 133714

write.table(merged_data_dedup_final, 
            "./table_STR_repeat_size_each_row_allele_EHv2.5.5_HTT_CAG_simplified_dedup_301219.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t"
            )
