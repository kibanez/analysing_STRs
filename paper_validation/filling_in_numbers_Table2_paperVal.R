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
table_diseases$year_of_birth = as.integer(table_diseases$year_of_birth)
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
  num_pid = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i]) %>% select(participant_id) %>% unique() %>% pull() %>% length()
  print(num_pid)
  
  # Age: median
  table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i]) %>% select(age) %>% pull() %>% mean() %>% print()
  
  # Age: 0-18 %
  age_0_18 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i], age > 0, age <= 18) %>% select(age) %>% pull() %>% length()) / num_pid
  
  # Age: 19-40%
  age_19_40 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i], age > 19, age <= 40) %>% select(age) %>% pull() %>% length()) / num_pid
  
  # Age: 41-60%
  age_41_60 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i], age > 41, age <= 60) %>% select(age) %>% pull() %>% length()) / num_pid
  
  # Age: 61-80%
  age_61_80 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i], age > 61, age <= 80) %>% select(age) %>% pull() %>% length()) / num_pid
  
  # Age: >80 %
  age_80 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main[i], age > 80) %>% select(age) %>% pull() %>% length()) / num_pid
  
  print(age_0_18)
  print(age_19_40)
  print(age_41_60)
  print(age_61_80)
  print(age_80)
}

# "Complex Parkinsonism" is different
table_diseases %>% filter(grepl("Complex Parkin", normalised_specific_disease)) %>% select(normalised_specific_disease) %>% unique() %>% print()
#normalised_specific_disease
# Complex Parkinsonism (includes pallido-pyramidal syndromes)
table_diseases %>% filter(grepl("Complex Parkin", normalised_specific_disease)) %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
# 139
table_diseases %>% filter(grepl("Complex Parkin", normalised_specific_disease)) %>% select(age) %>% pull() %>% mean() %>% print()

# Create the list of diseases in PILOT
l_diseases_pilot = c("Intellectual disability",
                     "Kabuki syndrome",
                     "Amyotrophic lateral sclerosis/motor neuron disease",
                     "Charcot-Marie-Tooth disease",
                     "Congenital muscular dystrophy",
                     "Congenital myopathy",
                     "Early onset dementia (encompassing fronto-temporal dementia and prion disease)",
                     "Early onset dystonia",
                     "Distal myopathies",
                     "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                     "Hereditary ataxia",
                     "Hereditary spastic paraplegia",
                     "Skeletal Muscle Channelopathies",
                     "Early onset and familial Parkinson's Disease")

for (i in 1:length(l_diseases_pilot)){
  # Check we are ONLY analysing the corresponding `disease`
  table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot[i]) %>% select(specificDisease) %>% unique() %>% print()
  
  # Number of UNIQUE participants
  num_pid = table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot[i]) %>% select(gelID) %>% unique() %>% pull() %>% length()
  print(num_pid)
  
  # Age: median
  table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot[i]) %>% select(age) %>% pull() %>% mean() %>% print()
  
  # Age: 0-18 %
  age_0_18 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot[i], age > 0, age <= 18) %>% select(age) %>% pull() %>% length()) / num_pid
  
  # Age: 19-40%
  age_19_40 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot[i], age > 19, age <= 40) %>% select(age) %>% pull() %>% length()) / num_pid
  
  # Age: 41-60%
  age_41_60 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot[i], age > 41, age <= 60) %>% select(age) %>% pull() %>% length()) / num_pid
  
  # Age: 61-80%
  age_61_80 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot[i], age > 61, age <= 80) %>% select(age) %>% pull() %>% length()) / num_pid
  
  # Age: >80 %
  age_80 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot[i], age > 80) %>% select(age) %>% pull() %>% length()) / num_pid
  
  print(age_0_18)
  print(age_19_40)
  print(age_41_60)
  print(age_61_80)
  print(age_80)
}


# Let's merge MAIN and PILOT
# This is the same loop
l_diseases_merge = intersect(l_diseases_main,
                             l_diseases_pilot)

for (i in 1:length(l_diseases_merge)){
  # Check we are ONLY analysing the corresponding `disease`
  table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i]) %>% select(normalised_specific_disease) %>% unique() %>% print()
  table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i]) %>% select(specificDisease) %>% unique() %>% print()
  
  # Number of UNIQUE participants
  num_pid_main = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i]) %>% select(participant_id) %>% unique() %>% pull() %>% length()
  num_pid_pilot = table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i]) %>% select(gelID) %>% unique() %>% pull() %>% length()
  total_num_pid = num_pid_main + num_pid_pilot 
  print(total_num_pid)
  
  # Age: median
  l_age_main = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i]) %>% select(age) %>% pull() 
  l_age_pilot = table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i]) %>% select(age) %>% pull()
  l_age_merged = c(l_age_main, l_age_pilot)
  print(mean(l_age_merged))
  
  # Age: 0-18 %
  main_age_0_18 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age >= 0, age <= 18) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_0_18 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i], age > 0, age <= 18) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
  merged_age_0_18 = sum(main_age_0_18, pilot_age_0_18) / total_num_pid
  
  # Age: 19-40%
  main_age_19_40 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age > 19, age <= 40) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_19_40 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i], age > 19, age <= 40) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_19_40 = sum(main_age_19_40, pilot_age_19_40) / total_num_pid
  
  # Age: 41-60%
  main_age_41_60 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age > 41, age <= 60) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_41_60 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i], age > 41, age <= 60) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_41_60 = sum(main_age_41_60, pilot_age_41_60) / total_num_pid
  
  # Age: 61-80%
  main_age_61_80 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age > 61, age <= 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_61_80 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i], age > 61, age <= 80) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_61_80 = sum(main_age_61_80, pilot_age_61_80) / total_num_pid
  
  # Age: >80 %
  main_age_80 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age > 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_80 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i], age > 80) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_80 = sum(main_age_80, pilot_age_80) / total_num_pid
  
  print(merged_age_0_18)
  print(merged_age_19_40)
  print(merged_age_41_60)
  print(merged_age_61_80)
  print(merged_age_80)
}


# This separately
extra_main = setdiff(l_diseases_main,
                     l_diseases_pilot)
# the only extra disease in MAIN not in PILOT is mitochondrial disorders
l_diseases_merge = c("Mitochondrial disorders")
for (i in 1:length(l_diseases_merge)){
  # Check we are ONLY analysing the corresponding `disease`
  table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i]) %>% select(normalised_specific_disease) %>% unique() %>% print()
  
  # Number of UNIQUE participants
  num_pid_main = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i]) %>% select(participant_id) %>% unique() %>% pull() %>% length()
  total_num_pid = num_pid_main 
  print(total_num_pid)
  
  # Age: median
  l_age_main = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i]) %>% select(age) %>% pull() 
  print(mean(l_age_main))
  
  # Age: 0-18 %
  main_age_0_18 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age >= 0, age <= 18) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_0_18 = sum(main_age_0_18) / total_num_pid
  
  # Age: 19-40%
  main_age_19_40 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age > 19, age <= 40) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_19_40 = sum(main_age_19_40) / total_num_pid
  
  # Age: 41-60%
  main_age_41_60 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age > 41, age <= 60) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_41_60 = sum(main_age_41_60) / total_num_pid
  
  # Age: 61-80%
  main_age_61_80 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age > 61, age <= 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_61_80 = sum(main_age_61_80) / total_num_pid
  
  # Age: >80 %
  main_age_80 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age > 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_80 = sum(main_age_80) / total_num_pid
  
  print(merged_age_0_18)
  print(merged_age_19_40)
  print(merged_age_41_60)
  print(merged_age_61_80)
  print(merged_age_80)
}



extra_pilot = setdiff(l_diseases_pilot,
                      l_diseases_main)


# Also Intellectual disability needs to be consider as `Kabuki` + ID
l_pid_ID_kabuki_main = table_diseases %>% 
  filter(normalised_specific_disease %in% c("Intellectual disability", "Kabuki syndrome")) %>% 
  select(participant_id) %>% 
  unique() %>%
  pull()

l_pid_ID_kabuki_pilot = table_diseases_pilot %>%
  filter(specificDisease %in% c("Intellectual disability", "Kabuki syndrome")) %>% 
  select(gelID) %>%
  unique() %>%
  pull()

l_pid_merged = c(l_pid_ID_kabuki_main,
                 l_pid_ID_kabuki_pilot)

length(l_pid_merged)
# 6731

# Mean age?
mean_age_main = table_diseases %>%
  filter(participant_id %in% l_pid_ID_kabuki_main) %>%
  select(age) %>%
  pull() %>%
  mean()

mean_age_pilot = table_diseases_pilot %>%
  filter(gelID %in% l_pid_ID_kabuki_pilot) %>%
  select(age) %>%
  pull() %>%
  mean()

mean_merged = mean(mean_age_main,
                   mean_age_pilot)
print(mean_merged)
# 13.33

# Age by distribution
main_age_0_18 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age >= 0, age <= 18) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_0_18 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age >= 0, age <= 18) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_0_18 = sum(main_age_0_18, pilot_age_0_18) / length(l_pid_merged)


main_age_19_40 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age >= 19, age <= 40) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_19_40 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age >= 19, age <= 40) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_19_40 = sum(main_age_19_40, pilot_age_19_40) / length(l_pid_merged)

main_age_41_60 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age >= 41, age <= 60) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_41_60 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age >= 41, age <= 60) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_41_60 = sum(main_age_41_60, pilot_age_41_60) / length(l_pid_merged)

main_age_61_80 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age >= 61, age <= 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_61_80 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age >= 61, age <= 80) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_61_80 = sum(main_age_61_80, pilot_age_61_80) / length(l_pid_merged)

main_age_more80 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age >= 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_more80 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age >= 80) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_more80 = sum(main_age_more80, pilot_age_more80) / length(l_pid_merged)
