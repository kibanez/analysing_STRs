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
# 12254  19

# Load pilot table
table_diseases_pilot = read.csv("table_diseases_enriched_PILOT_13diseases_enriched_popu.tsv",
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
dim(table_diseases_pilot)
# 660 13

# Removing platekey and genome assembly columns, since I'm interested in PID, and it's easier to deduplicate tables
table_diseases = table_diseases[-1]
table_diseases = table_diseases[-3]
table_diseases = unique(table_diseases)
dim(table_diseases)
# 12039  17

length(unique(table_diseases$participant_id))
# 11172 -- there are 867 duplicated pids, because they might have `disease_group` and `disease_subgroup`. 

# From here, I can see we need to focus only on the column names we want
table_diseases = table_diseases %>%
  select(participant_id, year_of_birth, normalised_specific_disease)
table_diseases_pilot = table_diseases_pilot %>%
  select(gelID, yearOfBirth, specificDisease)

# Define AGE, by using YOB
table_diseases$year_of_birth = as.integer(table_diseases$year_of_birth)
table_diseases = table_diseases %>%
  group_by(participant_id) %>%
  mutate(age = 2020 - year_of_birth) %>%
  ungroup() %>%
  as.data.frame()

table_diseases_pilot = table_diseases_pilot %>%
  group_by(gelID) %>%
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
                    "Kabuki syndrome",
                    "Complex Parkinsonism (includes pallido-pyramidal syndromes)")

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

# Let's merge MAIN and PILOT
# This is the same loop
l_diseases_merge = intersect(l_diseases_main,
                             l_diseases_pilot)

for (i in 1:length(l_diseases_merge)){
  # Check we are ONLY analysing the corresponding `disease`
  table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i]) %>% select(normalised_specific_disease) %>% unique() %>% print()
  table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i]) %>% select(specificDisease) %>% unique() %>% print()
  
  # Number of UNIQUE participants
  l_pid_disease_main = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i]) %>% select(participant_id) %>% unique() %>% pull()
  l_pid_disease_pilot = table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i]) %>% select(gelID) %>% unique() %>% pull()
  total_num_pid = length(l_pid_disease_main) + length(l_pid_disease_pilot)
  print(total_num_pid)
  
  # Age: median
  l_age_main = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i]) %>% select(age) %>% pull() 
  l_age_pilot = table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i]) %>% select(age) %>% pull()
  l_age_merged = c(l_age_main, l_age_pilot)
  print(mean(l_age_merged))
  
  # Age: 0-18 %
  main_age_0_18 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age %in% (0:18)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_0_18 = (table_diseases_pilot %>% filter(gelID %in% l_pid_disease_pilot & age %in% (0:18)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
  merged_age_0_18 = sum(main_age_0_18, pilot_age_0_18) / total_num_pid
  
  # Age: 19-40%
  main_age_19_40 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age %in% (19:40)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_19_40 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i], age %in% (19:40)) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_19_40 = sum(main_age_19_40, pilot_age_19_40) / total_num_pid
  
  # Age: 41-60%
  main_age_41_60 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age %in% (41:60)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_41_60 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i], age %in% (41:60)) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_41_60 = sum(main_age_41_60, pilot_age_41_60) / total_num_pid
  
  # Age: 61-80%
  main_age_61_80 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age %in% (61:80)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_61_80 = (table_diseases_pilot %>% filter(specificDisease %in% l_diseases_merge[i], age %in% (61:80)) %>% select(gelID) %>% unique() %>% pull() %>% length())
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
  main_age_0_18 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age %in% (0:18)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_0_18 = sum(main_age_0_18) / total_num_pid
  
  # Age: 19-40%
  main_age_19_40 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age %in% (19:40)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_19_40 = sum(main_age_19_40) / total_num_pid
  
  # Age: 41-60%
  main_age_41_60 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age %in% (41:60)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_41_60 = sum(main_age_41_60) / total_num_pid
  
  # Age: 61-80%
  main_age_61_80 = (table_diseases %>% filter(normalised_specific_disease %in% l_diseases_merge[i], age %in% (61:80)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
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

extra_main = setdiff(l_diseases_main,
                     l_diseases_pilot)

extra_pilot = setdiff(l_diseases_pilot,
                      l_diseases_main)
extra_main = extra_main[-4]

for (i in 1:length(extra_main)){
  # Check we are ONLY analysing the corresponding `disease`
  table_diseases %>% filter(normalised_specific_disease %in% extra_main[i]) %>% select(normalised_specific_disease) %>% unique() %>% print()
  table_diseases_pilot %>% filter(specificDisease %in% extra_pilot[i]) %>% select(specificDisease) %>% unique() %>% print()
  
  # Number of UNIQUE participants
  num_pid_main = table_diseases %>% filter(normalised_specific_disease %in% extra_main[i]) %>% select(participant_id) %>% unique() %>% pull() %>% length()
  num_pid_pilot = table_diseases_pilot %>% filter(specificDisease %in% extra_pilot[i]) %>% select(gelID) %>% unique() %>% pull() %>% length()
  total_num_pid = num_pid_main + num_pid_pilot 
  print(total_num_pid)
  
  # Age: median
  l_age_main = table_diseases %>% filter(normalised_specific_disease %in% extra_main[i]) %>% select(age) %>% pull() 
  l_age_pilot = table_diseases_pilot %>% filter(specificDisease %in% extra_pilot[i]) %>% select(age) %>% pull()
  l_age_merged = c(l_age_main, l_age_pilot)
  print(mean(l_age_merged))
  
  # Age: 0-18 %
  main_age_0_18 = (table_diseases %>% filter(normalised_specific_disease %in% extra_main[i], age %in% (0:18)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_0_18 = (table_diseases_pilot %>% filter(specificDisease %in% extra_pilot[i], age %in% (0:18)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
  merged_age_0_18 = sum(main_age_0_18, pilot_age_0_18) / total_num_pid
  
  # Age: 19-40%
  main_age_19_40 = (table_diseases %>% filter(normalised_specific_disease %in% extra_main[i], age %in% (19:40)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_19_40 = (table_diseases_pilot %>% filter(specificDisease %in% extra_pilot[i], age %in% (19:40)) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_19_40 = sum(main_age_19_40, pilot_age_19_40) / total_num_pid
  
  # Age: 41-60%
  main_age_41_60 = (table_diseases %>% filter(normalised_specific_disease %in% extra_main[i], age %in% (41:60)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_41_60 = (table_diseases_pilot %>% filter(specificDisease %in% extra_pilot[i], age %in% (41:60)) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_41_60 = sum(main_age_41_60, pilot_age_41_60) / total_num_pid
  
  # Age: 61-80%
  main_age_61_80 = (table_diseases %>% filter(normalised_specific_disease %in% extra_main[i], age %in% (61:80)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_61_80 = (table_diseases_pilot %>% filter(specificDisease %in% extra_pilot[i], age %in% (61:80)) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_61_80 = sum(main_age_61_80, pilot_age_61_80) / total_num_pid
  
  # Age: >80 %
  main_age_80 = (table_diseases %>% filter(normalised_specific_disease %in% extra_main[i], age > 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  pilot_age_80 = (table_diseases_pilot %>% filter(specificDisease %in% extra_pilot[i], age > 80) %>% select(gelID) %>% unique() %>% pull() %>% length())
  merged_age_80 = sum(main_age_80, pilot_age_80) / total_num_pid
  
  print(merged_age_0_18)
  print(merged_age_19_40)
  print(merged_age_41_60)
  print(merged_age_61_80)
  print(merged_age_80)
}

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
main_age_0_18 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age %in% c(0:18)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_0_18 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age %in% c(0:18)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_0_18 = sum(main_age_0_18, pilot_age_0_18) / length(l_pid_merged)

main_age_19_40 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age %in% c(19:40)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_19_40 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age %in% c(19:40)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_19_40 = sum(main_age_19_40, pilot_age_19_40) / length(l_pid_merged)

main_age_41_60 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age %in% c(41:60)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_41_60 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age %in% c(41:60)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_41_60 = sum(main_age_41_60, pilot_age_41_60) / length(l_pid_merged)

main_age_61_80 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age %in% c(61:80)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_61_80 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age %in% c(61:80)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_61_80 = sum(main_age_61_80, pilot_age_61_80) / length(l_pid_merged)

main_age_more80 = (table_diseases %>% filter(participant_id %in% l_pid_ID_kabuki_main, age >= 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_more80 = (table_diseases_pilot %>% filter(gelID %in% l_pid_ID_kabuki_pilot, age >= 80) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_more80 = sum(main_age_more80, pilot_age_more80) / length(l_pid_merged)


# total of participants
# Load main table
table_diseases = read.csv("table_diseases_enriched_popu_includingSkeletalMuscleChan.tsv",
                          stringsAsFactors = F, 
                          header = T,
                          sep = "\t")
dim(table_diseases)
# 12254  19

# Load pilot table
table_diseases_pilot = read.csv("table_diseases_enriched_PILOT_13diseases_enriched_popu.tsv",
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
dim(table_diseases_pilot)
# 660 13

# males vs females
# From here, I will select the correpsonding columns I want to work with
table_diseases = table_diseases %>%
  select(participant_id, participant_phenotypic_sex, year_of_birth, participant_ethnic_category)
table_diseases_pilot = table_diseases_pilot %>%
  select(gelID, sex, yearOfBirth, specificDisease)

table_diseases = unique(table_diseases)
table_diseases_pilot = unique(table_diseases_pilot)

# Define AGE, by using YOB
table_diseases$year_of_birth = as.integer(table_diseases$year_of_birth)
table_diseases = table_diseases %>%
  group_by(participant_id) %>%
  mutate(age = 2020 - year_of_birth) %>%
  ungroup() %>%
  as.data.frame()

table_diseases_pilot = table_diseases_pilot %>%
  group_by(gelID) %>%
  mutate(age = 2020 - yearOfBirth) %>%
  ungroup() %>%
  as.data.frame()

# How many male vs female?
# Main
l_main_female = table_diseases %>% filter(participant_phenotypic_sex %in% "Female") %>% select(participant_id) %>% unique() %>% pull() 
length(l_main_female)
# 4749
l_main_male = table_diseases %>% filter(participant_phenotypic_sex %in% "Male") %>% select(participant_id) %>% unique() %>% pull() 
length(l_main_male)
# 6423

# Pilot
l_pilot_female = table_diseases_pilot %>% filter(sex %in% "female") %>% select(gelID) %>% unique() %>% pull()
length(l_pilot_female)
# 297
l_pilot_male = table_diseases_pilot %>% filter(sex %in% "male") %>% select(gelID) %>% unique() %>% pull()
length(l_pilot_male)
# 348


l_merged_female = c(l_main_female, l_pilot_female)
l_merged_male = c(l_main_male, l_pilot_male)
length(l_merged_female)
# 5046
length(l_merged_male)
# 6771

# Age mean and breakdown
# Mean age?
mean_age_main_female = table_diseases %>%
  filter(participant_id %in% l_main_female) %>%
  select(age) %>%
  pull() %>%
  mean()

mean_age_main_male = table_diseases %>%
  filter(participant_id %in% l_main_male) %>%
  select(age) %>%
  pull() %>%
  mean()


mean_age_pilot_female = table_diseases_pilot %>%
  filter(gelID %in% l_pilot_female) %>%
  select(age) %>%
  pull() %>%
  mean()

mean_age_pilot_male = table_diseases_pilot %>%
  filter(gelID %in% l_pilot_male) %>%
  select(age) %>%
  pull() %>%
  mean()


mean_merged_female = mean(mean_age_main_female,
                   mean_age_pilot_female)

mean_merged_male = mean(mean_age_main_male,
                          mean_age_pilot_male)

print(mean_merged_female)
# 29.02
print(mean_merged_male)
# 25.9

# for females
main_age_0_18 = (table_diseases %>% filter(participant_id %in% l_main_female, age %in% c(0:18)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_0_18 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_female, age %in% c(0:18)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_0_18 = sum(main_age_0_18, pilot_age_0_18) / length(l_merged_female)

main_age_19_40 = (table_diseases %>% filter(participant_id %in% l_main_female, age %in% c(19:40)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_19_40 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_female, age %in% c(19:40)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_19_40 = sum(main_age_19_40, pilot_age_19_40) / length(l_merged_female)

main_age_41_60 = (table_diseases %>% filter(participant_id %in% l_main_female, age %in% c(41:60)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_41_60 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_female, age %in% c(41:60)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_41_60 = sum(main_age_41_60, pilot_age_41_60) / length(l_merged_female)

main_age_61_80 = (table_diseases %>% filter(participant_id %in% l_main_female, age %in% c(61:80)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_61_80 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_female, age %in% c(61:80)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_61_80 = sum(main_age_61_80, pilot_age_61_80) / length(l_merged_female)

main_age_more80 = (table_diseases %>% filter(participant_id %in% l_main_female, age >= 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_more80 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_female, age >= 80) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_more80 = sum(main_age_more80, pilot_age_more80) / length(l_merged_female)

# for males
main_age_0_18 = (table_diseases %>% filter(participant_id %in% l_main_male, age %in% c(0:18)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_0_18 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_male, age %in% c(0:18)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_0_18 = sum(main_age_0_18, pilot_age_0_18) / length(l_merged_male)

main_age_19_40 = (table_diseases %>% filter(participant_id %in% l_main_male, age %in% c(19:40)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_19_40 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_male, age %in% c(19:40)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_19_40 = sum(main_age_19_40, pilot_age_19_40) / length(l_merged_male)

main_age_41_60 = (table_diseases %>% filter(participant_id %in% l_main_male, age %in% c(41:60)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_41_60 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_male, age %in% c(41:60)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_41_60 = sum(main_age_41_60, pilot_age_41_60) / length(l_merged_male)

main_age_61_80 = (table_diseases %>% filter(participant_id %in% l_main_male, age %in% c(61:80)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_61_80 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_male, age %in% c(61:80)) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_61_80 = sum(main_age_61_80, pilot_age_61_80) / length(l_merged_male)

main_age_more80 = (table_diseases %>% filter(participant_id %in% l_main_male, age >= 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
pilot_age_more80 = (table_diseases_pilot %>% filter(gelID %in% l_pilot_male, age >= 80) %>% select(gelID) %>% unique() %>% pull() %>% length()) 
merged_age_more80 = sum(main_age_more80, pilot_age_more80) / length(l_merged_male)

# Let's compute now the ethnicity reported (just for Main, since Pilot is lacking that info)
l_ethnias = unique(table_diseases$participant_ethnic_category)
length(l_ethnias)
# 18

# Let's see mean and breakdown age across existing 18 diff ethnias
for (i in 1:length(l_ethnias)){
  # Check we are ONLY analysing the corresponding `disease`
  table_diseases %>% filter(participant_ethnic_category %in% l_ethnias[i]) %>% select(participant_ethnic_category) %>% unique() %>% print()
  
  # Number of UNIQUE participants
  num_pid_main = table_diseases %>% filter(participant_ethnic_category %in% l_ethnias[i]) %>% select(participant_id) %>% unique() %>% pull() %>% length()
  total_num_pid = num_pid_main 
  print(total_num_pid)
  
  # Age: median
  l_age_main = table_diseases %>% filter(participant_ethnic_category %in% l_ethnias[i]) %>% select(age) %>% pull() 
  print(mean(l_age_main))
  
  # Age: 0-18 %
  main_age_0_18 = (table_diseases %>% filter(participant_ethnic_category %in% l_ethnias[i], age %in% (0:18)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_0_18 = sum(main_age_0_18) / total_num_pid
  
  # Age: 19-40%
  main_age_19_40 = (table_diseases %>% filter(participant_ethnic_category %in% l_ethnias[i], age %in% (19:40)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_19_40 = sum(main_age_19_40) / total_num_pid
  
  # Age: 41-60%
  main_age_41_60 = (table_diseases %>% filter(participant_ethnic_category %in% l_ethnias[i], age %in% (41:60)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_41_60 = sum(main_age_41_60) / total_num_pid
  
  # Age: 61-80%
  main_age_61_80 = (table_diseases %>% filter(participant_ethnic_category %in% l_ethnias[i], age %in% (61:80)) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_61_80 = sum(main_age_61_80) / total_num_pid
  
  # Age: >80 %
  main_age_80 = (table_diseases %>% filter(participant_ethnic_category %in% l_ethnias[i], age > 80) %>% select(participant_id) %>% unique() %>% pull() %>% length())
  merged_age_80 = sum(main_age_80) / total_num_pid
  
  print(merged_age_0_18)
  print(merged_age_19_40)
  print(merged_age_41_60)
  print(merged_age_61_80)
  print(merged_age_80)
  
}