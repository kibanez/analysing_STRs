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

# The age cut-off if only for : ALS, charcot, early onset dementia, early onset dystonia, complex parkinson, hereditary spastic para-legaa and easy; onset familial parkin
table_diseases_select = table_diseases %>%
  filter(normalised_specific_disease %in% l_diseases_main)
table_diseases_select = rbind(table_diseases_select,
                              table_diseases %>% filter(grepl("Complex Parkin", normalised_specific_disease)))

table_diseases_pilot_select = table_diseases_pilot %>%
  filter(specificDisease %in% l_diseases_pilot)

# We need to take ONLY ADULTS for all diseases for except for 
# Keep everyone for `hereditary ataxia`, and `Intelectual.. and kabuki`




# Load latest clinical data for these diseases for the paper - April 2021
clin_data = read.csv("./table_diseases_for_table2_Main_and_Pilot_14785_PIDs_all_adults_and_paediatrics_withouth_Mitochondrial.tsv",
                     sep = "\t",
                     stringsAsFactors = F)
dim(clin_data)
# 28300 26

# Define AGE, by using YOB - keeping year of the analysis
clin_data = clin_data %>%
  group_by(participant_id) %>%
  mutate(age = 2020 - year_of_birth) %>%
  ungroup() %>%
  as.data.frame()

# Let's recode the ethnicity, simplifying it
table_diseases$participant_ethnic_category = recode(table_diseases$participant_ethnic_category,
                                             "White: British"= "White",
                                             "White White: British"= "White",
                                             "White: Any other White background"="White",
                                             "Asian or Asian British: Pakistani"="Asian",
                                             "Asian or Asian British: Indian"="Asian",
                                             "Asian or Asian British: Any other Asian background"="Asian",
                                             "Black or Black British: African"="Black",
                                             "Other Ethnic Groups: Any other ethnic group"="Other", 
                                             "Mixed: Any other mixed background"="Mixed",
                                             "White: Irish"="White",
                                             "Asian or Asian British: Bangladeshi"="Asian",
                                             "Mixed: White and Asian"="Mixed",
                                             "Mixed: White and Black Caribbean"="Mixed",
                                             "Black or Black British: Caribbean"="Black",
                                             "Mixed: White and Black African"="Mixed",
                                             "Black or Black British: Any other Black background"="Black",
                                             "Other Ethnic Groups: Chinese"="Other")
# Defining NA's as `Not stated`
which_na = which(is.na(table_diseases$participant_ethnic_category))
table_diseases$participant_ethnic_category[which_na] = "Not Stated"

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

# We need to keep the PIDs that have these diseases (already selected)
# All participants (adults and paediatrics) except ONLY ADULTS for:
# `Hereditary spastic paraplegia`
# `Early onset and familial Parkinson''s Disease`
# `Complex Parkin`
# `Early onset dystonia`
# `Early onset dementia`
# `Amyotrophic lateral sclerosis`
# `Charcot-Marie-Tooth disease`

# First, retrieve ONLY adults for these diseases
clin_data_adults_only = clin_data %>%
  filter((grepl("Hereditary spastic paraplegia", diseases_list, ignore.case = T) & age >= 18) |
         (grepl("Early onset and familial Parkinson''s Disease", diseases_list, ignore.case = T) & age >= 18) |
         (grepl("Complex Parkin", diseases_list, ignore.case = T) & age >= 18) |
         (grepl("Early onset dystonia", diseases_list, ignore.case = T) & age >= 18) |
         (grepl("Early onset dementia", diseases_list, ignore.case = T) & age >= 18) |
         (grepl("Amyotrophic lateral sclerosis", diseases_list, ignore.case = T) & age >= 18) |
         (grepl("Charcot-Marie-Tooth disease", diseases_list, ignore.case = T) & age >= 18))

# Check average age
summary(clin_data_adults_only$age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18.00   35.00   50.00   50.02   65.00   99.00 

# First, remove from `clin_data` the following PIDs having at least any disease for ADULTS ONLY
clin_data = clin_data %>%
  filter(!grepl("Hereditary spastic paraplegia", diseases_list, ignore.case = T))
clin_data = clin_data %>%
  filter(!grepl("Early onset and familial Parkinson''s Disease", diseases_list, ignore.case = T))
clin_data = clin_data %>%
  filter(!grepl("Complex Parkin", diseases_list, ignore.case = T))
clin_data = clin_data %>%
  filter(!grepl("Early onset dystonia", diseases_list, ignore.case = T))
clin_data = clin_data %>%
  filter(!grepl("Early onset dementia", diseases_list, ignore.case = T))
clin_data = clin_data %>%
  filter(!grepl("Amyotrophic lateral sclerosis", diseases_list, ignore.case = T))
clin_data = clin_data %>%
  filter(!grepl("Charcot-Marie-Tooth disease", diseases_list, ignore.case = T))

# Merge ADULTS only for the specified disease-list + the rest
clin_data = rbind(clin_data,
                  clin_data_adults_only)  
clin_data = unique(clin_data)
dim(clin_data)
# 27680  27

length(unique(clin_data$participant_id))
# 14108
length(unique(clin_data$latest_platekey))
# 14108

# Select randomly 13331 from above list
set.seed(512)
l_PID = sample(unique(clin_data$participant_id), 13331)

clin_data = clin_data %>%
  filter(participant_id %in% l_PID)

# Mean age
summary(clin_data$age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.0    11.0    18.0    27.8    43.0   102.0 

# Gender
aver = clin_data %>% select(participant_id, participant_phenotypic_sex) %>% unique()
table(aver$participant_phenotypic_sex)
#female        Female Indeterminate          male          Male 
#269          5534             1           320          7207 
# 5803 females
# 7527 males

l_independent_diseases_adults = c("Hereditary spastic paraplegia",
                                  "'Early onset and familial Parkinson''s Disease'",
                                  "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                                  "Early onset dystonia",
                                  "Early onset dementia",
                                  "Amyotrophic lateral sclerosis or motor neuron disease",
                                  "Charcot-Marie-Tooth disease")

l_independent_diseases_pilot_adults = c("Hereditary spastic paraplegia",
                                        "Early onset and familial Parkinson's Disease",
                                        "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                                        "Early onset dystonia",
                                        "Early onset dementia (encompassing fronto-temporal dementia and prion disease)",
                                        "Amyotrophic lateral sclerosis/motor neuron disease",
                                        "Charcot-Marie-Tooth disease")



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

# From this list of 13K pids, I'll save it
write.table(l_pid_all, "./list_13331_PIDs_table2.txt", row.names = F, col.names = F, quote = F)

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
aver = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main) %>% select(participant_phenotypic_sex, participant_id)  %>% unique() 
table(aver$participant_phenotypic_sex)
#Female   Male 
#5477   7071

aver2 = table_diseases %>% filter(grepl(l_diseases_main_extra, normalised_specific_disease)) %>% select(participant_phenotypic_sex, participant_id) %>% unique() 
table(aver2$participant_phenotypic_sex)
#Female   Male 
#62     77 

aver_pilot = table_diseases_pilot %>% filter(specificDisease %in% l_diseases_pilot) %>% select(sex, gelID) %>% unique() 
table(aver_pilot$sex)
#female   male 
#297    348 

# Ethnicity
l_main_eth = table_diseases %>% filter(normalised_specific_disease %in% l_diseases_main) %>% select(participant_ethnic_category)  %>% pull() 
l_main_eth_extra = table_diseases %>% filter(grepl(l_diseases_main_extra, normalised_specific_disease)) %>% select(participant_ethnic_category) %>% pull() 
l_eth_all = c(l_main_eth,
              l_main_eth_extra)
table(l_eth_all)
#Asian      Black      Mixed Not Stated      Other      White 
#1235        264        384       2109        174       9132 


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
                     
# Only adults for
# Hereditary spastic paraplegia, 
# Early onset and familial Parkinson's Disease
# Complex Parkinsonism (includes pallido-pyramidal syndromes)
# Early onset dystonia, 
# Early onset dementia, 
# Amyotrophic lateral sclerosis or motor neuron disease
# Charcot-Marie-Tooth disease
l_independent_diseases_adults = c("Hereditary spastic paraplegia",
                           "'Early onset and familial Parkinson''s Disease'",
                           "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                           "Early onset dystonia",
                           "Early onset dementia",
                           "Amyotrophic lateral sclerosis or motor neuron disease",
                           "Charcot-Marie-Tooth disease")


l_independent_diseases_pilot_adults = c("Hereditary spastic paraplegia",
                                 "Early onset and familial Parkinson's Disease",
                                 "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                                 "Early onset dystonia",
                                 "Early onset dementia (encompassing fronto-temporal dementia and prion disease)",
                                 "Amyotrophic lateral sclerosis/motor neuron disease",
                                 "Charcot-Marie-Tooth disease")

for (i in 1:length(l_independent_diseases_adults)){
  print(l_independent_diseases_adults[i])
  # Unique PIDs
  l_main_pid = table_diseases %>% filter(normalised_specific_disease %in% l_independent_diseases_adults[i], adult.paediatric %in% "Adult") %>% select(participant_id) %>% unique() %>% pull() 
  l_pilot_pid = table_diseases_pilot %>% filter(specificDisease %in% l_independent_diseases_pilot_adults[i], adult.paediatric %in% "Adult") %>% select(gelID) %>% unique() %>% pull() 
  
  l_pid_all = unique(c(l_main_pid,
                       l_pilot_pid))
  print(length(l_pid_all))
  
  # Age of all
  l_main_age = table_diseases %>% filter(normalised_specific_disease %in% l_independent_diseases_adults[i], adult.paediatric %in% "Adult") %>% select(age) %>% pull() 
  l_pilot_age = table_diseases_pilot %>% filter(specificDisease %in% l_independent_diseases_pilot_adults[i], adult.paediatric %in% "Adult") %>% select(age) %>% pull() 
  
  l_age_all = c(l_main_age,
                l_pilot_age)
  print(summary(l_age_all))
  
  # Gender
  l_main_gender = table_diseases %>% filter(normalised_specific_disease %in% l_independent_diseases_adults[i], adult.paediatric %in% "Adult") %>% select(participant_phenotypic_sex)  %>% pull() 
  l_pilot_gender = table_diseases_pilot %>% filter(specificDisease %in% l_independent_diseases_pilot_adults[i], adult.paediatric %in% "Adult") %>% select(sex) %>% pull() 
  
  l_gender_all = c(l_main_gender,
                   l_pilot_gender)
  
  print(table(l_gender_all))
  
  # Ethnicity
  l_main_eth = table_diseases %>% filter(normalised_specific_disease %in% l_independent_diseases_adults[i], adult.paediatric %in% "Adult") %>% select(participant_ethnic_category)  %>% pull() 
  print(table(l_main_eth))
  
  
  
}

# From `clin_data` which is already filtered out by 13331 and only adults for some specific disease
l_independent_diseases = c("Hereditary ataxia",
                           "Hereditary spastic paraplegia",
                           "Early onset and familial Parkinson's Disease",
                           "Complex Parkin",
                           "Early onset dystonia",
                           "Early onset dementia",
                           "Amyotrophic lateral sclerosis",
                           "Charcot-Marie-Tooth disease",
                           "Ultra-rare undescribed monogenic disorders",
                           "Intellectual disability",
                           "Kabuki syndrome",
                           "Congenital myopathy",
                           "Distal myopathies",
                           "Congenital muscular dystrophy",
                           "Skeletal Muscle Channelopathies")
            
for (i in 1:length(l_independent_diseases)){
  print(l_independent_diseases[i])
  
  l_disease_PID = clin_data %>% filter(grepl(l_independent_diseases[i], diseases_list)) %>% select(participant_id) %>% unique() %>% pull()
  print(length(l_disease_PID))
  
  # Age of all
  clin_data %>% filter(participant_id %in% l_disease_PID) %>% select(age) %>% summary() %>% print()
  
  # Gender
  aver = clin_data %>% filter(participant_id %in% l_disease_PID) %>% select(participant_id, participant_phenotypic_sex) %>% unique()
  print(table(aver$participant_phenotypic_sex))
}

