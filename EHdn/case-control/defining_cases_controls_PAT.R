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
# 154  10

length(unique(pilot_cases$plateKey))
# 153


main_cases = main_clin_data %>%
  filter(year_of_birth < 2000,
         participant_type %in% "Proband",
         !grepl("Neurology", main_clin_data$disease_group),
         genome_build %in% "GRCh37")
dim(main_cases)
# 16711  28

length(unique(main_cases$platekey))
# 1410


