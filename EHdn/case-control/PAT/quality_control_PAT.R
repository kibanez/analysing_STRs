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
  filter(plateKey %in% l_pilot_cases) %>%
  select(plateKey, sex, yearOfBirth)
dim(selected_pilot_cases)
# 194  3

selected_main_cases = main_cases %>%
  filter(platekey %in% l_main_cases) %>%
  select(platekey, participant_phenotypic_sex, year_of_birth)
dim(selected_main_cases)
# 59432  3

selected_main_cases = unique(selected_main_cases)
dim(selected_main_cases)
# 1427  3

colnames(selected_pilot_cases) = c("platekey", "sex", "age")
colnames(selected_main_cases) = c("platekey", "sex", "age")

selected_cases = rbind(selected_pilot_cases,
                       selected_main_cases)

dim(selected_cases)
# 1621 3

# Pilot and Main controls
selected_pilot_controls = pilot_controls %>%
  filter(plateKey %in% l_pilot_controls) %>%
  select(plateKey, sex, yearOfBirth)
dim(selected_pilot_controls)
# 1145  3

selected_pilot_controls = unique(selected_pilot_controls)
dim(selected_pilot_controls)
# 1113  3

selected_main_controls = main_controls %>%
  filter(platekey %in% l_main_controls) %>%
  select(platekey, participant_phenotypic_sex, year_of_birth)
dim(selected_main_controls)
# 48648  3

selected_main_controls = unique(selected_main_controls)
dim(selected_main_controls)
# 2371  3

colnames(selected_pilot_controls) = c("platekey", "sex", "age")
colnames(selected_main_controls) = c("platekey", "sex", "age")

selected_controls = rbind(selected_pilot_controls,
                          selected_main_controls)
dim(selected_controls)
# 3483  3

# Define each `selected_cases` and `selected_controls` with the group name
selected_cases$group = rep("case", length(selected_cases$platekey))
selected_controls$group = rep("control", length(selected_controls$platekey))

merged_pat = rbind(selected_cases,
                   selected_controls)

dim(merged_pat)
# 109419 4


# Plot age distribution
ggplot()


ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot()
