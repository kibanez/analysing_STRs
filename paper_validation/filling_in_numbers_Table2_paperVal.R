# Objective: R script to fill out Table2 of the paper

# libraries
library(dplyr)
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1

# Set directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

# Load main table
table_diseases = read.csv("table_diseases_enriched_popu_includingSkeletalMuscleChan.tsv",
                          stringsAsFactors = F, 
                          header = T,
                          sep = "\t")
dim(table_diseases)
# 11801  16

# Load pilot table
table_diseases_pilot = read.csv("table_diseases_enriched_PILOT_13diseases_enriched_popu.tsv",
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
dim(table_diseases_pilot)
# 660 13

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

# Counting number of participants, age (mean and distribution)

# Create the list of diseases in MAIN
l_diseases_main = c("Intellectual disability",
                    "Amyotrophic lateral sclerosis or motor neuron disease", 
                    "Charcot-Marie-Tooth disease", 
                    "Congenital muscular dystrophy",
                    "Congenital myopathy", 
                    "Early onset dementia", 
                    "Early onset dystonia", 
                    "Distal myopathies", 
                    "Hereditary ataxia", 
                    "Hereditary spastic paraplegia", 
                    "Skeletal Muscle Channelopathies",
                    "'Early onset and familial Parkinson''s Disease'",
                    "Mitochondrial disorders",
                    "Kabuki syndrome")

for (i in 1:length(l_diseases_main)){
  # Check we are ONLY analysing the corresponding `disease`
  table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i]) %>% select(normalised_specific_disease) %>% unique() %>% print()
  
  # Number of UNIQUE participants
  table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i]) %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
  
  # Age: median
  table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i]) %>% select(age) %>% unique() %>% pull() %>% length() %>% print()
}

# "Complex Parkinsonism" is different
table_diseases %>% filter(grepl("Complex Parkin", normalised_specific_disease)) %>% select(normalised_specific_disease) %>% unique() %>% print()
#normalised_specific_disease
# Complex Parkinsonism (includes pallido-pyramidal syndromes)
table_diseases %>% filter(grepl("Complex Parkin", normalised_specific_disease)) %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
# 139


table_diseases_pilot %>% filter()