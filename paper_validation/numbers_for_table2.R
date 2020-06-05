# Objective: count the total number of participants, %gender. %ethnicity 
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
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

# Load table with the diagnostics 
# Main table
table_diseases = read.csv("table_diseases_enriched_popu_includingSkeletalMuscleChan_and_ultra-rare.tsv",
                          stringsAsFactors = F, 
                          header = T,
                          sep = "\t")
dim(table_diseases)
# 13868  19

# Pilot table
table_diseases_pilot = read.csv("table_diseases_enriched_PILOT_13diseases_enriched_popu.tsv",
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
dim(table_diseases_pilot)
# 660  13

# Define AGE, by using YOB
table_diseases = table_diseases %>%
  group_by(participant_id) %>%
  mutate(age = 2020 - year_of_birth) %>%
  ungroup() %>%
  as.data.frame()

table_diseases_pilot = table_diseases_pilot %>%
  group_by(plateKey) %>%
  mutate(age = 2020 - yearOfBirth) %>%
  ungroup() %>%
  as.data.frame()

# Define adult or paediatric
table_diseases = mutate(table_diseases, adult.paediatric = ifelse(age < 18, "Paediatric", "Adult"))
table_diseases_pilot = mutate(table_diseases_pilot, adult.paediatric = ifelse(age < 18, "Paediatric", "Adult"))
