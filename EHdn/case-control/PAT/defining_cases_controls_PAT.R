# Objective: define and track here the modus operandi when selecting samples for cases and controls for PILOT adult ataxia
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/case-control/analysis/PAT/")

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

# I've seen that the `phenotyping_table` file if not well formatted, and for each `specific_disease` there are disease_group and disease_subgroup == blank
# we need to clean them up
phenotyping_table = phenotyping_table %>% 
  filter(!disease_group %in% "")

dim(phenotyping_table)
# 2599  3

# enrich pilot_clin_data with associations with disease_group and disease_subgroup
pilot_clin_data = inner_join(pilot_clin_data,
                            phenotyping_table,
                            by = c("specificDisease" = "specific_disease"))
dim(pilot_clin_data)
# 146591  12

pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 2105  12

main_clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T)
dim(main_clin_data)
# 1124633  28

# Cases
# - ONLY probands
# - YOB < 2000 (i.e. adults)
# - specific disease in ataxia
# - for `main` take GRCh37

pilot_cases = pilot_clin_data %>%
  filter(yearOfBirth < 2000, 
         grepl("[Aa]taxia", specificDisease), 
         biological_relation_to_proband %in% "Proband")
dim(pilot_cases)
# 117  12

length(unique(pilot_cases$plateKey))
# 117

l_pilot_cases = unique(pilot_cases$plateKey)

pilot_controls = pilot_clin_data %>%
  filter(biological_relation_to_proband %in% "Proband",
         !grepl("[Nn][Ee][Uu][Rr][Oo]",disease_group),
         !plateKey %in% l_pilot_cases,
         yearOfBirth < 2000)
dim(pilot_controls)
# 892  12

length(unique(pilot_controls$plateKey))
# 865

# Controls
# - ONLY probands
# - rare disease germline genomes (exclude for now cancer genomes)
# - `disease_group` NOT in neurology
# - for `main` take GRCh37

main_cases = main_clin_data %>%
  filter(participant_type %in% "Proband",
         year_of_birth < 2000,
         genome_build %in% "GRCh37",
         grepl("[Aa]taxia", specific_disease))
dim(main_cases)
# 4350  28

length(unique(main_cases$plate_key))
# 97

l_main_cases = unique(main_cases$plate_key)

main_controls = main_clin_data %>%
  filter(participant_type %in% "Proband",
         programme %in% "Rare Diseases",
         genome_build %in% "GRCh37",
         year_of_birth < 2000,
         !plate_key %in% l_main_cases,
         !grepl("[Nn][Ee][Uu][Rr][Oo]", main_clin_data$disease_group))
dim(main_controls)
# 16693  28

length(unique(main_controls$plate_key))
# 1391

# Writing individual files
l_pilot_cases = unique(pilot_cases$plateKey)
write.table(l_pilot_cases, "input/pilot_117_cases.txt", quote = F, row.names = F, col.names = F)

l_pilot_controls = unique(pilot_controls$plateKey)
write.table(l_pilot_controls, "input/pilot_865_controls.txt", quote = F, row.names = F, col.names = F)

l_main_cases = unique(main_cases$plate_key)
write.table(l_main_cases, "input/main_97_cases.txt", quote = F, row.names = F, col.names = F)
# which(l_main_cases %in% l_platekey_both)
write.table(l_main_cases[-which(l_main_cases %in% l_platekey_both)], "input/main_79_cases.txt", quote = F, row.names = F, col.names = F)


l_main_controls = unique(main_controls$plate_key)
write.table(l_main_controls, "input/main_1391_controls.txt", quote = F, row.names = F, col.names = F)
# which(l_main_controls %in% l_platekey_both)
write.table(l_main_controls[-which(l_main_controls %in% l_platekey_both)], "input/main_1199_controls.txt", quote = F, row.names = F, col.names = F)


# Merged CASE and CONTROL files
#l_cases = unique(c(l_pilot_cases,
#                   l_main_cases))

l_cases = unique(c(l_pilot_cases,
                   l_main_cases[-which(l_main_cases %in% l_platekey_both)]))

#l_controls = unique(c(l_pilot_controls,
#                      l_main_controls))

l_controls = unique(c(l_pilot_controls,
                      l_main_controls[-which(l_main_controls %in% l_platekey_both)]))

#write.table(l_cases, "input/merged_pilot_main_214_cases.txt", quote = F, row.names = F, col.names = F)
#write.table(l_controls, "input/merged_pilot_main_2256_controls.txt", quote = F, row.names = F, col.names = F)

# only b37
write.table(l_cases, "input/merged_pilot_main_196_cases.txt", quote = F, row.names = F, col.names = F)
write.table(l_controls, "input/merged_pilot_main_2064_controls.txt", quote = F, row.names = F, col.names = F)


# Let's create now the `manifest` file
# We need to merge all STR profiles for all case-control samples together into a multi-sample STR profile. 
# This requires us to create the manifest file manifest.tsv describing the dataset. 
# The manifest file contains columns for sample identifier, case-control status, and path to the associated STR profile for each sample.

#sample1 case    str-profiles/sample1.str_profile.json
#sample2 case    str-profiles/sample2.str_profile.json
#sample3 case    str-profiles/sample3.str_profile.json
#sample4 control str-profiles/sample4.str_profile.json
#sample5 control str-profiles/sample5.str_profile.json
#sample6 control str-profiles/sample6.str_profile.json
#sample7 control str-profiles/sample7.str_profile.json

cases_df = data.frame(platekey = l_cases,
                      group = rep("case", length(l_cases)))

cases_df = cases_df %>%
  group_by(platekey) %>%
  mutate(path = paste(paste("/genomes/scratch/kgarikano/GEL_STR/EHdn/case-control/analysis/PAT/str-profiles", platekey, sep = "/"), "_EHdeNovo.str_profile.json", sep = "")) %>%
  ungroup() %>%
  as.data.frame()

controls_df = data.frame(platekey = l_controls,
                      group = rep("control", length(l_controls)))

controls_df = controls_df %>%
  group_by(platekey) %>%
  mutate(path = paste(paste("/genomes/scratch/kgarikano/GEL_STR/EHdn/case-control/analysis/PAT/str-profiles", platekey, sep = "/"), "_EHdeNovo.str_profile.json", sep = "")) %>%
  ungroup() %>%
  as.data.frame()

merged_df = rbind(cases_df,
                  controls_df)

dim(merged_df)
# 2470  3
# 2260  3

# QC check - there should not be duplicated platekeys
length(merged_df$platekey)
# 2470
# 2260

write.table(merged_df,
            "input/manifest_PAT_only_b37.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

save.image("PAT_case_control_environment.Rdata")

# run quality control checks
source("~/git/analysing_STRs/EHdn/case-control/functions/quality_control.R")
plotting_age_distribution(environment_file = "~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/case-control/analysis/PAT/PAT_case_control_environment.Rdata", 
                          working_directory = "~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/case-control/analysis/PAT/")





