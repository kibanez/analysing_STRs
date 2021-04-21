# Objective: create the 4 tables for the diagnostic results section
# April 2021
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
table_diseases = read.csv("april2021_tableDiseases_12686_individuals_main.tsv",
                          stringsAsFactors = F, 
                          header = T,
                          sep = "\t")
dim(table_diseases)
# 13298  16

# Pilot table
table_diseases_pilot = read.csv("april2021_tableDiseases_645_individuals_pilot.tsv",
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
dim(table_diseases_pilot)
# 660  11

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

# Let's compute the denominator of number of individuals included in Panel A






# Let's define now the 4 subtables for the purpose of the paper

# TABLE A. ONLY INCLUDING ADULTS (I.E. >= 18 IN 2020, EXCEPT FXN WHERE WE INCLUDE CHILDREN), USING FULL-MUTATION CUTOFF THRESHOLD  
# (OR YOU CAN PRODUCE A TABLE USING THE PREMUTATION CUTOFF, BUT I SUSPOECT IT WILL BE VERY NOISY AND WOULD NOT REFELCT THE THRESHOLDS THAT PANELAPP IS CURRENTLY USING)

# Change on 22nd April: 
# add ADULT patients recruited under: 'Ultra-rare undescribed monogenic disorders' (SPECIFIC DISEASE), 
# only the ones that have been applied one or more of the following PANELS:  
# Amyotrophic lateral sclerosis/motor neuron disease, 
# Hereditary neuropathy, 
# Early onset dementia (encompassing fronto-temporal dementia and prion disease), 
# Parkinson Disease and Complex Parkinsonism, 
# Early onset dystonia, 
# Hereditary spastic paraplegia, 
# Hereditary ataxia

# select diseases we are interested for TABLE A
table_a = table_diseases %>%
  filter(normalised_specific_disease %in% c("Amyotrophic lateral sclerosis or motor neuron disease", 
                                            "Charcot-Marie-Tooth disease",
                                            "Early onset dementia", 
                                            "Early onset dystonia", 
                                            "Complex Parkinsonism", 
                                            "Hereditary ataxia", 
                                            "Hereditary spastic paraplegia",
                                            "'Early onset and familial Parkinson''s Disease'"))
dim(table_a)
# 3518  21

# Complex parkinsonism is missing here
table_a = rbind(table_a,
                table_diseases %>%
                  filter(grepl("[Cc]omplex [Pp]arkin", table_diseases$normalised_specific_disease)))
dim(table_a)
# 3659  21

table_HA =  table_a %>%
  filter(normalised_specific_disease %in% "Hereditary ataxia")

table_a = table_a %>%
  filter(!normalised_specific_disease%in% "Hereditary ataxia",
         adult.paediatric %in% "Adult")

# Let's define list of diseases for Table A, as we have done for the genes
l_diseases_tableA = unique(table_a$normalised_specific_disease)
length(l_diseases_tableA)
# 8

# Let's select panels from `Ultra-rare disorder``
# List panels
list_panels_part2 = c("Amyotrophic lateral sclerosis/motor neuron disease",
                      " Amyotrophic lateral sclerosis/motor neuron disease",
                      "Hereditary neuropathy",
                      " Hereditary neuropathy",
                      "Early onset dementia (encompassing fronto-temporal dementia and prion disease)",
                      " Early onset dementia (encompassing fronto-temporal dementia and prion disease)",
                      "Parkinson Disease and Complex Parkinsonism",
                      " Parkinson Disease and Complex Parkinsonism",
                      "Early onset dystonia",
                      " Early onset dystonia",
                      "Hereditary spastic paraplegia",
                      " Hereditary spastic paraplegia",
                      "Hereditary ataxia",
                      " Hereditary ataxia")

# let's make life simple: split into rows panel info
table_panels_row = table_diseases %>% 
  select(plate_key.x, participant_id, normalised_specific_disease, disease_sub_group, disease_group, panel_list, adult.paediatric) %>%
  mutate(panels = strsplit(as.character(panel_list), ",")) %>%
  unnest(panels) %>%
  as.data.frame()
table_panels_row = unique(table_panels_row)
dim(table_panels_row)
# 58449  8
table_panels_row$participant_id = as.character(table_panels_row$participant_id)

# Function that checks if any of the items in list of characters 1 does exist in list of characters 2
any_exist <- function(list1, list2) {
  list2_sep = unique(strsplit(list2, ',')[[1]])
  for (i in list1){
    if (i %in% list2_sep){
      return(TRUE)
    }
  }
  return(FALSE)
}


# select diseases we are interested for TABLE A - part2
table_a_part2 = table_diseases %>%
  filter(normalised_specific_disease %in% "Ultra-rare undescribed monogenic disorders",
         adult.paediatric %in% "Adult")
dim(table_a_part2)
# 725  21

# Let's take only those that have any of the panels specified in part2 within list_panels
table_a_part2 = table_a_part2 %>%
  group_by(participant_id) %>%
  mutate(any_panel_in_listpanels_part2 = any_exist(list_panels_part2, panel_list)) %>%
  ungroup() %>%
  as.data.frame()
dim(table_a_part2)
# 725  22 

# just take the ones including any of the panels suggested
table_a_part2 = table_a_part2 %>%
  filter(any_panel_in_listpanels_part2)
dim(table_a_part2)
# 49  22

table_a_part2 = table_a_part2[,-19]

table_a = rbind(table_a,
                table_HA,
                table_a_part2)
table_a = unique(table_a)
dim(table_a)
# 3398  18

# How many PIDs?
length(unique(table_a$participant_id))
# 3279
# How many platekeys?
length(unique(table_a$plate_key.x))
# 3279

# select diseases we are interested for TABLE A - PILOT
table_a_pilot = table_diseases_pilot %>%
  filter(specificDisease %in% c("Amyotrophic lateral sclerosis/motor neuron disease",
                                "Charcot-Marie-Tooth disease",
                                "Early onset dementia (encompassing fronto-temporal dementia and prion disease)",
                                "Early onset dystonia",
                                "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                                "Hereditary ataxia",
                                "Hereditary spastic paraplegia",
                                "Early onset and familial Parkinson's Disease"))
dim(table_a_pilot)
# 418  15

table_a_pilot_HA = table_a_pilot %>%
  filter(specificDisease %in% "Hereditary ataxia")

table_a_pilot = table_a_pilot %>%
  filter(!specificDisease %in% "Hereditary ataxia",
         adult.paediatric %in% "Adult")

table_a_pilot = rbind(table_a_pilot,
                      table_a_pilot_HA)
table_a_pilot = unique(table_a_pilot)
dim(table_a_pilot)
# 411 13

# How many PIDs are in the Pilot?
length(unique(table_a_pilot$plateKey))
# 401
length(unique(table_a_pilot$gelID))
# 401

# Let's select the interesting columns for Table 2
table_a = table_a %>% select(participant_id, plate_key.x, rare_diseases_family_id, participant_phenotypic_sex, year_of_birth, normalised_specific_disease)
table_a_pilot = table_a_pilot %>% select(gelID, plateKey, gelFamilyId.x, sex, yearOfBirth, specificDisease)
colnames(table_a_pilot) = colnames(table_a)

# Let's merge MAIN and PILOT
panel_a = rbind(table_a,
                table_a_pilot)
panel_a = unique(panel_a)
dim(panel_a)
# 3697  6

length(unique(panel_a$participant_id))
# 3680
panel_a$panel = rep("A", length(panel_a$participant_id))

################################################################################################################################################################
# TABLE B
# We need to take the list of PIDs from `list_2459_PIDs_ID_and_others_as_panels.txt`

# load list of 2449 PIDs 
l_complex_ID_group2 = read.table("./list_2743_PIDs_ID_and_others_as_panels.txt", stringsAsFactors = F)
l_complex_ID_group2 = l_complex_ID_group2$V1
length(l_complex_ID_group2)
# 2743

table_b = table_diseases %>%
  filter(participant_id %in% l_complex_ID_group2)
dim(table_b)
# 3132  19

length(unique(table_b$plate_key.x))
# 2743

length(unique(table_b$participant_id))
# 2743

panel_b = table_b %>% select(participant_id, plate_key.x, rare_diseases_family_id, participant_phenotypic_sex, year_of_birth, normalised_specific_disease)
panel_b$panel = rep("B", length(panel_b$participant_id))
panel_b = unique(panel_b)
dim(panel_b)
# 2851  7

length(unique(panel_b$participant_id))
# 2743
################################################################################################################################################################
# TABLE C
# patients presenting with intellectual disability and or a neuromuscular phenotype were analysed for DMPK

# MAIN
table_c = table_diseases %>%
  filter(normalised_specific_disease %in% c("Intellectual disability",
                                            "Kabuki syndrome",
                                            "Congenital muscular dystrophy",
                                            "Congenital myopathy",
                                            "Skeletal Muscle Channelopathies",
                                            "Distal myopathies"))
dim(table_c)
# 7695  19

# PILOT
table_c_pilot = table_diseases_pilot %>%
  filter(specificDisease %in% c("Intellectual disability",
                                "Kabuki syndrome",
                                "Congenital muscular dystrophy",
                                "Congenital myopathy",
                                "Skeletal Muscle Channelopathies",
                                "Distal myopathies"))
dim(table_c_pilot)
# 242  13

table_c = table_c %>% select(participant_id, plate_key.x, rare_diseases_family_id, participant_phenotypic_sex, year_of_birth, normalised_specific_disease)
table_c_pilot = table_c_pilot %>% select(gelID, plateKey, gelFamilyId.x, sex, yearOfBirth, specificDisease)
colnames(table_c_pilot) = colnames(table_c)

panel_c = rbind(table_c,
                table_c_pilot)
panel_c = unique(panel_c)
dim(panel_c)
# 7592  6

panel_c$panel = rep("C", length(panel_c$participant_id))

# How many PIDs
length(unique(panel_c$participant_id))
# 7586
length(unique(panel_c$plate_key.x))
# 7586

################################################################################################################################################################
# TABLE D. only including children recruited under ID (using >55. as cutoff)	
# FMR1
# Intellectual disability	

# MAIN
table_d = table_diseases %>%
  filter(normalised_specific_disease %in% c("Intellectual disability",
                                            "Kabuki syndrome"))
dim(table_d)
# 6890  19

# PILOT
table_d_pilot = table_diseases_pilot %>%
  filter(specificDisease %in% c("Intellectual disability",
                                "Kabuki syndrome"))
dim(table_d_pilot)
# 161  13

table_d = table_d %>% select(participant_id, plate_key.x, rare_diseases_family_id, participant_phenotypic_sex, year_of_birth, normalised_specific_disease)
table_d_pilot = table_d_pilot %>% select(gelID, plateKey, gelFamilyId.x, sex, yearOfBirth, specificDisease)
colnames(table_d_pilot) = colnames(table_d)

panel_d = rbind(table_d,
                table_d_pilot)
panel_d = unique(panel_d)
dim(panel_d)
# 6731  6

# PIDs?
length(unique(panel_d$participant_id))
# 6731

