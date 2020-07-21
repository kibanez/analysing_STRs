# Objective: count how many probands we do have
# Count how many probands MINUS neuro we do have
# For EHv3.2.2
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2019-02-29)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/pileup_100Kg/")

# All loci tables should have the same number of pids...probands, etc.
# ACHTUNG! case-controls have all unique PIDs, but they might be related
l_unrelated_main = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/60k_HWE_30k_random_unrelated_participants.txt",
                              stringsAsFactors = F,
                              header = F)
l_unrelated_main = l_unrelated_main$V1
length(l_unrelated_main)
# 38344

df_pilot = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/data_freeze_Pilot_LK_RESEARCH/pedigree.csv",
                      stringsAsFactors = F,
                      header = T,
                      sep = ",")
dim(df_pilot)
# 17258  34

l_pid_unrelated_pilot = df_pilot %>%
  filter(isProband %in% "TRUE") %>%
  select(gelId) %>%
  unique() %>%
  pull()
length(l_pid_unrelated_pilot)
# 2279

df_pilot2 = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/b37_genomes_2019-12-11_11-07-20.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(df_pilot2)
# 4828  5

l_unrelated_pilot = df_pilot2 %>%
  filter(participant_id %in% l_pid_unrelated_pilot) %>%
  select(plate_key) %>%
  unique() %>%
  pull()
length(l_unrelated_pilot)
# 2278

# Let's write into a file the list of unrelated MAIN and PILOT platekeys

l_unrelated = c(l_unrelated_main,
                l_unrelated_pilot)
length(l_unrelated)
# 40622

write.table(l_unrelated, 
            "~/Documents/STRs/clinical_data/l_unrelated_40622_platekeys_MAIN_and_PILOT.txt",
            quote = F,
            row.names = F,
            col.names = F)

# From here, with the list of unrelated MAIN and PILOT platekeys, let's compute the total number of probands, probands-neuro
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_research_cohort/clin_data_merged_V5:V9.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(clin_data)
# 2061403  31 

# Probands ONLY
probands_unrelated = clin_data %>%
  filter(plate_key %in% l_unrelated, participant_type %in% "Proband") %>%
  select(plate_key) %>%
  unique() %>%
  pull()
length(probands_unrelated)
# 13076

# there are still some relatives...let's remove them
which_relative = clin_data %>%
  filter(plate_key %in% probands_unrelated, participant_type %in% "Relative") %>%
  select(plate_key) %>%
  unique() %>%
  pull()

index_to_remove_as_relative = which(probands_unrelated %in% which_relative)
probands_unrelated = probands_unrelated[-index_to_remove_as_relative]
length(probands_unrelated)
# 13036

# Probands ONLY - (minus or except) NEUROLOGY
# Remove from here participants that are not part of NEURO
probands_not_neuro_unrelated = clin_data %>%
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group), plate_key %in% l_unrelated, participant_type %in% "Proband") %>%
  select(plate_key) %>%
  unique() %>%
  pull()
length(probands_not_neuro_unrelated)
# 8299

# Again, there are platekeys that are still neuro... and proband...let's filter them out
# Ths is because there might be participants that have been recruited as Father/Mother, and eventually has been assigned as Affected 
which_relative = clin_data %>%
  filter(plate_key %in% probands_not_neuro_unrelated, participant_type %in% "Relative") %>%
  select(plate_key) %>%
  unique() %>%
  pull()

index_to_remove_as_relative = which(probands_not_neuro_unrelated %in% which_relative)
probands_not_neuro_unrelated = probands_not_neuro_unrelated[-index_to_remove_as_relative]
length(probands_not_neuro_unrelated)
# 8270

which_neuro = clin_data %>%
  filter(plate_key %in% probands_not_neuro_unrelated, grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group)) %>%
  select(plate_key) %>%
  unique() %>%
  pull()

index_to_remove_as_neuro = which(probands_not_neuro_unrelated %in% which_neuro)
probands_not_neuro_unrelated = probands_not_neuro_unrelated[-index_to_remove_as_neuro]
length(probands_not_neuro_unrelated)
# 8198

# Write them into files
write.table(probands_unrelated,
            "list_probands_unrelated_13036_platekeys.txt",
            quote = F,
            row.names = F,
            col.names = F)

write.table(probands_not_neuro_unrelated,
            "list_probands_unrelated_NOT_NEURO_8198_platekeys.txt",
            quote = F,
            row.names = F,
            col.names = F)
