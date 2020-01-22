# Objective: define and track here the modus operandi when selecting samples for cases and controls for PILOT adult ataxia
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/case-control/")

# in this case, the first analysis we will perform is called PAT: pilot adult ataxia 
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           sep = "\t",
                           stringsAsFactors = F,
                           header = T)
dim(pilot_clin_data)
# 4974  10


# the table with `disease_group` and `specific_disease` associations
phenotyping_table = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/phenotyping_v140_2019-09-13_15-26-02.tsv",
                             sep = "\t",
                             stringsAsFactors = F,
                             header = T)
dim(phenotyping_table)
# 2632  3

# enrich pilot_clin_data with associations with disease_group and disease_subgroup
pilot_clin_data = left_join(pilot_clin_data,
                            phenotyping_table,
                            by = c("specificDisease" = "specific_disease"))
dim(pilot_clin_data)
# 150912  12

pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 5587  12

main_clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T)
dim(main_clin_data)
# 1124633  28

# Cases
# - ONLY probands
# - YOB < 2000 (i.e. adults)
# - specific disease in ataxia OR neuropathies
# - for `main` take GRCh37

pilot_cases = pilot_clin_data %>%
  filter(yearOfBirth < 2000, 
         (grepl("neuropathies", specificDisease) | grepl("ataxia", specificDisease)), 
         biological_relation_to_proband %in% "Proband")
dim(pilot_cases)
# 154  12

length(unique(pilot_cases$plateKey))
# 153

pilot_controls = pilot_clin_data %>%
  filter(biological_relation_to_proband %in% "Proband",
         !grepl("[Nn][Ee][Uu][Rr][Oo]",disease_group))
dim(pilot_controls)
# 2274  12

length(unique(pilot_controls$plateKey))
# 1957

# Controls
# - ONLY probands
# - rare disease germline genomes (exclude for now cancer genomes)
# - `disease_group` NOT in neurology
# - for `main` take GRCh37

main_cases = main_clin_data %>%
  filter(participant_type %in% "Proband",
         year_of_birth < 2000,
         (grepl("neuropathies", specific_disease) | grepl("ataxia", specific_disease)))
dim(main_cases)
# 45779  28

length(unique(main_cases$platekey))
# 1014

main_controls = main_clin_data %>%
  filter(participant_type %in% "Proband",
         programme %in% "Rare Diseases",
         genome_build %in% "GRCh37",
         !grepl("[Nn][Ee][Uu][Rr][Oo]", main_clin_data$disease_group))
dim(main_controls)
# 48648  28

length(unique(main_controls$platekey))
# 2371


# Writing individual files
l_pilot_cases = unique(pilot_cases$plateKey)
write.table(l_pilot_cases, "./PAT/pilot_153_cases.txt", quote = F, row.names = F, col.names = F)

l_pilot_controls = unique(pilot_controls$plateKey)
write.table(l_pilot_controls, "./PAT/pilot_1957_controls.txt", quote = F, row.names = F, col.names = F)

l_main_cases = unique(main_cases$platekey)
write.table(l_main_cases, "./PAT/main_1014_cases.txt", quote = F, row.names = F, col.names = F)

l_main_controls = unique(main_controls$platekey)
write.table(l_main_controls, "./PAT/main_2371_controls.txt", quote = F, row.names = F, col.names = F)

# Merged CASE and CONTROL files
l_cases = unique(c(l_pilot_cases,
                   l_main_cases))

l_controls = unique(c(l_pilot_controls,
                      l_main_controls))

write.table(l_cases, "./PAT/merged_pilot_main_1167_cases.txt", quote = F, row.names = F, col.names = F)
write.table(l_controls, "./PAT/merged_pilot_main_4328_controls.txt", quote = F, row.names = F, col.names = F)




