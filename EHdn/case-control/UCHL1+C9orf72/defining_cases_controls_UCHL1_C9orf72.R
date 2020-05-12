# Objective: define and track here the modus operandi when selecting samples for cases and controls for PILOT adult ataxia
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.9.0/case-control/analysis/UCHL1_C9orf72/")

# Load list of cases
l_cases = read.table("./input/list_4_genomes_UCLH1.txt", stringsAsFactors = F)
l_cases = l_cases$V1


# Load Main data
main_clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_300320.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T)
dim(main_clin_data)
# 1124633  31

# Cases - `l_cases`
l_main_cases = l_cases

# Controls
# - ONLY probands
# - rare disease germline genomes (exclude for now cancer genomes)
# - `disease_group` NOT in neurology
# - for `main` take GRCh38
main_controls = main_clin_data %>%
  filter(participant_type %in% "Proband",
         programme %in% "Rare Diseases",
         genome_build %in% "GRCh38",
         year_of_birth < 2000,
         !plate_key %in% l_main_cases,
         !grepl("[Nn][Ee][Uu][Rr][Oo]", main_clin_data$disease_group))
dim(main_controls)
# 124282  31

length(unique(main_controls$plate_key))
# 11355

# Writing into files
write.table(l_main_cases, "input/main_4_cases.txt", quote = F, row.names = F, col.names = F)

l_main_controls = unique(main_controls$plate_key)
write.table(l_main_controls, "input/main_11355_controls.txt", quote = F, row.names = F, col.names = F)

# Let's create now the `manifest` file
cases_df = data.frame(platekey = l_main_cases,
                      group = rep("case", length(l_cases)))

cases_df = cases_df %>%
  group_by(platekey) %>%
  mutate(path = paste(paste("/home/kgarikano/GEL_STR/EHdn/case-control/EHdn_v0.9.0/analysis/UCHL1_C9orf72/str-profiles", platekey, sep = "/"), "_EHdeNovo.str_profile.json", sep = "")) %>%
  ungroup() %>%
  as.data.frame()

controls_df = data.frame(platekey = l_main_controls,
                         group = rep("control", length(l_main_controls)))

controls_df = controls_df %>%
  group_by(platekey) %>%
  mutate(path = paste(paste("/home/kgarikano/GEL_STR/EHdn/case-control/EHdn_v0.9.0/analysis/UCHL1_C9orf72/str-profiles", platekey, sep = "/"), "_EHdeNovo.str_profile.json", sep = "")) %>%
  ungroup() %>%
  as.data.frame()

merged_df = rbind(cases_df,
                  controls_df)

dim(merged_df)
# 11359  3

# QC check - there should not be duplicated platekeys
length(merged_df$platekey)
# 11359

write.table(merged_df,
            "input/manifest_UCHL1.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)


# There are some genomes that have not been included in EHdn calculation
real_list = read.table("../test_ALS/list_11226_cases.txt", stringsAsFactors = F)
real_list = real_list$V1
# 11226

merged_df = merged_df %>%
  filter(group %in% "case" | (platekey %in% real_list & group %in% "control"))
dim(merged_df)
# 11225  3

write.table(merged_df,
            "input/manifest_UCHL1_real.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)


save.image("UCHL1_C9orf72_case_control_environment.Rdata")

# run quality control checks
source("~/git/analysing_STRs/EHdn/case-control/functions/quality_control.R")
plotting_age_distribution(environment_file = "~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.9.0/case-control/analysis/test_ALS/test_ALS_case_control_environment.Rdata", 
                          working_directory = "~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.9.0/case-control/analysis/test_ALS/")
