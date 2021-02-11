# Objective: select genomes sequenced on GRCh37 from MAIN and PILOT programmes
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)
library(reshape)
library(tidyverse)

setwd("/Users/kibanez/Documents/STRs/ANALYSIS/EHdn/GRCh37/input/")

clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V10_and_Pilot_programmes.tsv",
                     stringsAsFactors = F, 
                     header = T,
                     sep = "\t")
dim(clin_data)
# 2101385  24

# CASES selection - Main GRCh37 + Pilot
# `normalised_specific_disease` in ataxia OR hereditary spastic paraplegia
cases_pilot_ha = clin_data %>%
  filter(programme %in% "RD Pilot", grepl("Hereditary ataxia", diseases_list))

cases_pilot_hsp = clin_data %>%
  filter(programme %in% "RD Pilot", grepl("Hereditary spastic paraplegia", diseases_list))

cases_pilot = unique(rbind(cases_pilot_ha,
                           cases_pilot_hsp))
dim(cases_pilot)
# 243  24

cases_main_ha = clin_data %>%
  filter(programme %in% "Rare Diseases", grepl("Hereditary ataxia", diseases_list), genome_build %in% "GRCh37")

cases_main_hsp = clin_data %>%
  filter(programme %in% "Rare Diseases", grepl("Hereditary spastic paraplegia", diseases_list), genome_build %in% "GRCh37")

cases_main = unique(rbind(cases_main_ha,
                          cases_main_hsp))
dim(cases_main)
# 405  24

merged_cases = unique(rbind(cases_pilot,
                            cases_main))
dim(merged_cases)
# 648  24

l_platekeys_cases = unique(merged_cases$platekey)
length(l_platekeys_cases)
# 471

# CONTROL selection - Main GRCh37 + Pilot
#normalised_specific_disease as N/A (no disease)
#disease_group not in Neurology/NEUROLOGY/Neurodevelopmental diseases
#check that the genomes considered as cases are not within this group
#for main dataset - select only GRCh37 genomes
controls_pilot = clin_data %>%
  filter(programme %in% "RD Pilot", !grepl("[Nn][Ee][Uu][Rr][Oo]", diseasegroup_list))

controls_main = clin_data %>%
  filter(programme %in% "Rare Diseases", !grepl("[Nn][Ee][Uu][Rr][Oo]", diseasegroup_list), !platekey %in% l_platekeys_cases, genome_build %in% "GRCh37")

merged_controls = unique(rbind(controls_pilot,
                               controls_main))
dim(merged_controls)
# 18173  24

l_platekeys_controls = unique(merged_controls$platekey)
length(l_platekeys_controls)
# 12374