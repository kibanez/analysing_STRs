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


