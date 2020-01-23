# Objective: to see age distribution across case and control datasets
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/case-control/PAT/")

# Load PAT case-control environment
load("PAT_case_control_environment.Rdata")

# Pilot and Main cases - taking platekey, sex and age
selected_pilot_cases = pilot_cases %>% 
  select(plateKey, sex, yearOfBirth)

selected_main_cases = main_cases %>%
  select(platekey, participant_phenotypic_sex, year_of_birth)

colnames(selected_pilot_cases) = c("platekey", "sex", "age")
colnames(selected_main_cases) = c("platekey", "sex", "age")

selected_cases = rbind(selected_pilot_cases,
                       selected_main_cases)

dim(selected_cases)
# 59626  3

# Pilot and Main controls
selected_pilot_controls = pilot_controls %>%
  select(plateKey, sex, yearOfBirth)

selected_main_controls = main_controls %>%
  select(platekey, participant_phenotypic_sex, year_of_birth)

colnames(selected_pilot_controls) = c("platekey", "sex", "age")
colnames(selected_main_controls) = c("platekey", "sex", "age")

selected_controls = rbind(selected_pilot_controls,
                          selected_main_controls)
