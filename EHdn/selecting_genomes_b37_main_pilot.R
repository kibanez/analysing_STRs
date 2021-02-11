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

# CONTROL selection - Main GRCh37 + Pilot