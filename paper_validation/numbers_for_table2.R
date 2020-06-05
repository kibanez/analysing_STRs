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

# List of normalised specific diseases
l_diseases_main = c("Intellectual disability",
               "Amyotrophic lateral sclerosis or motor neuron disease", 
               "Charcot-Marie-Tooth disease", "Congenital muscular dystrophy",
               "Congenital myopathy", "Early onset dementia", "Early onset dystonia", 
               "Distal myopathies", "Complex Parkinsonism", "Hereditary ataxia", 
               "Hereditary spastic paraplegia", "Skeletal Muscle Channelopathies",
               "'Early onset and familial Parkinson''s Disease'",
               "Mitochondrial disorders",
               "Kabuki syndrome",
               "Ultra-rare undescribed monogenic disorders")

l_diseases_main_extra = "Complex Parkin"

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

# How many unique PIDs are in MAIN and PILOT cohorts with these diseases?
l_main_pid = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main) %>% select(participant_id) %>% unique() %>% pull() 
# 12548
l_main_pid_extra = table_diseases %>% filter(grepl(l_diseases_main_extra, normalised_specific_disease)) %>% select(participant_id) %>% unique() %>% pull() 
# 139
l_pilot_pid = table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot) %>% select(gelID) %>% unique() %>% pull() 
# 645

l_pid_all = unique(c(l_main_pid,
                l_main_pid_extra,
                l_pilot_pid))
length(l_pid_all)
# 13331 

# Age of all
l_main_age = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main) %>% select(age) %>% pull() 
l_main_age_extra = table_diseases %>% filter(grepl(l_diseases_main_extra, normalised_specific_disease)) %>% select(age)  %>% pull() 
l_pilot_age = table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot) %>% select(age) %>% pull() 

l_age_all = c(l_main_age,
              l_main_age_extra,
              l_pilot_age)
summary(l_age_all)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.00    9.00   17.00   27.34   45.00  101.00 

# Gender
l_main_gender = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main) %>% select(participant_phenotypic_sex)  %>% pull() 
l_main_gender_extra = table_diseases %>% filter(grepl(l_diseases_main_extra, normalised_specific_disease)) %>% select(participant_phenotypic_sex) %>% pull() 
l_pilot_gender = table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot) %>% select(sex) %>% pull() 

l_gender_all = c(l_main_gender,
                 l_main_gender_extra,
                 l_pilot_gender)

table(l_gender_all)
# female Female   male   Male 
# 302   5817    358   7481 

# female =  6119
# male = 7839

# total = 13958

# Let's calculate across independent diseases
l_independent_diseases = c("Hereditary ataxia",
                           "Hereditary spastic paraplegia",
                           "'Early onset and familial Parkinson''s Disease'",
                           "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                           "Early onset dystonia",
                           "Early onset dementia",
                           "Amyotrophic lateral sclerosis or motor neuron disease",
                           "Charcot-Marie-Tooth disease",
                           "Ultra-rare undescribed monogenic disorders",
                           "Intellectual disability",
                           "Kabuki syndrome",
                           "Congenital myopathy",
                           "Distal myopathies",
                           "Congenital muscular dystrophy",
                           "Skeletal Muscle Channelopathies")


l_independent_diseases_pilot = c("Hereditary ataxia",
                     "Hereditary spastic paraplegia",
                     "Early onset and familial Parkinson's Disease",
                     "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                     "Early onset dystonia",
                     "Early onset dementia (encompassing fronto-temporal dementia and prion disease)",
                     "Amyotrophic lateral sclerosis/motor neuron disease",
                     "Charcot-Marie-Tooth disease",
                     "",
                     "Intellectual disability",
                     "Kabuki syndrome",
                     "Congenital myopathy",
                     "Distal myopathies",
                     "Congenital muscular dystrophy",
                     "Skeletal Muscle Channelopathies")
                     
            
for (i in 1:length(l_independent_diseases)){
  print(l_independent_diseases[i])
  # Unique PIDs
  l_main_pid = table_diseases %>% filter(normalised_specific_disease %in% l_independent_diseases[i]) %>% select(participant_id) %>% unique() %>% pull() 
  l_pilot_pid = table_diseases_pilot %>% filter(specificDisease %in% l_independent_diseases_pilot[i]) %>% select(gelID) %>% unique() %>% pull() 
  
  l_pid_all = unique(c(l_main_pid,
                       l_pilot_pid))
  print(length(l_pid_all))
  
  # Age of all
  l_main_age = table_diseases %>% filter(normalised_specific_disease %in% l_independent_diseases[i]) %>% select(age) %>% pull() 
  l_pilot_age = table_diseases_pilot %>% filter(specificDisease %in% l_independent_diseases_pilot[i]) %>% select(age) %>% pull() 
  
  l_age_all = c(l_main_age,
                l_pilot_age)
  print(summary(l_age_all))

  # Gender
  l_main_gender = table_diseases %>% filter(normalised_specific_disease %in% l_independent_diseases[i]) %>% select(participant_phenotypic_sex)  %>% pull() 
  l_pilot_gender = table_diseases_pilot %>% filter(specificDisease %in% l_independent_diseases_pilot[i]) %>% select(sex) %>% pull() 
  
  l_gender_all = c(l_main_gender,
                   l_pilot_gender)
  
  print(table(l_gender_all))
  
  
}

