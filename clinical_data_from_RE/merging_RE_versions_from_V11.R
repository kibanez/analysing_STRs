# Objective: have a unique file which merges all RE versions (Vxx) with the most important info for RD
# PID, FAMID, Platekey, RE Version
# So, whenever we want to know when the sample was removed or is included, we can easily know this
# This is a simplified version from RE V11 (to avoid the loading of all versions)
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
setwd("~/Documents/STRs/clinical_data/clinical_data/")

# Load the catalog-file with MAIN RE V1 to V11 + PILOT
clin_data = read.csv("./merged_RE_releases_and_Pilot_RD_and_Cancer_PID_FID_platekey.tsv",
                     stringsAsFactors = F,
                     sep = "\t",
                     header = T)
dim(clin_data)
# 920109  6

clin_data = unique(clin_data)
dim(clin_data)
# 175060 6

clin_release_V12 = read.csv("./rd_genomes_all_data_240521_V12.tsv",
                            stringsAsFactors = F,
                            sep = "\t",
                            header = T)
dim(clin_release_V12)
# 2105636 36

clin_data = clin_data %>%
  group_by(participant_id) %>%
  mutate(list_re_versions2 = ifelse(participant_id %in% clin_release_V12$participant_id, paste(list_re_versions, "RE_V12", sep = ","), list_re_versions)) %>%
  ungroup() %>%
  as.data.frame()

clin_data = unique(clin_data)
dim(clin_data)
# 175060  7

# Include extra genomes in V12 not in previous versions (83)
l_pid_new_in_V12 = setdiff(clin_release_V12$participant_id, clin_data$participant_id)
to_include = clin_release_V12 %>%
  filter(participant_id %in% l_pid_new_in_V12) %>%
  select(participant_id, platekey, rare_diseases_family_id, type, genome_build)

to_include$list_re_versions2 = rep("RE_V12", length(to_include$participant_id))

# remove previous `list_re_versions`
clin_data = clin_data[,-5]

clin_data = rbind(clin_data,
                  to_include)

# Write into a file
write.table(clin_data,
            "~/Documents/STRs/clinical_data/clinical_data/merged_RE_releases_and_Pilot_RD_and_Cancer_PID_FID_platekey_up_to_RE_V12.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
