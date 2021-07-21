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

# Load  family_history table
is_fami = read.csv("./table_family_history_for_11k_PIDs_MAIN_PROGRAMME.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(is_fami)
# 67006  11

is_fami_pilot = read.csv("./table_family_history_for_11k_PIDs_PILOT_PROGRAMME.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(is_fami_pilot)
# 1308  11

is_fami_pilot = is_fami_pilot %>%
  group_by(gelFamilyId.x) %>%
  mutate(is_fam_pilot = ifelse(num_affected_members > 1, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame() %>%
  unique()

# Load table with the diagnostics 
# Main table
table_diseases = read.csv("./table_diseases_enriched_popu_includingSkeletalMuscleChan_and_ultra-rare.tsv",
                          stringsAsFactors = F, 
                          header = T,
                          sep = "\t")
dim(table_diseases)
# 13868  19

# Pilot table -  22nd April 2021: we include now ultra-rare corresponding diseases from Pilot
table_diseases_pilot = read.csv("table_diseases_enriched_PILOT_13diseases_22April2021.tsv",
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
dim(table_diseases_pilot)
# 831  11

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

table_a_part2 = table_a_part2[,-22]

table_a = rbind(table_a,
                table_HA,
                table_a_part2)
table_a = unique(table_a)
dim(table_a)
# 3398  21

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
# 418  13

table_a_pilot_HA = table_a_pilot %>%
  filter(specificDisease %in% "Hereditary ataxia")

table_a_pilot = table_a_pilot %>%
  filter(!specificDisease %in% "Hereditary ataxia",
         adult.paediatric %in% "Adult")

table_a_pilot_ultra = table_diseases_pilot %>%
  filter(specificDisease %in% c("Unknown disorder",
                                "All recognised syndromes and those with suggestive features"),
         adult.paediatric %in% "Adult")
dim(table_a_pilot_ultra)
# 133 13

# let's make life simple: split into rows panel info
table_panels_row_pilot = table_diseases_pilot %>% 
  select(plateKey, gelID, specificDisease, panel_list, adult.paediatric) %>%
  mutate(panels = strsplit(as.character(panel_list), ",")) %>%
  unnest(panels) %>%
  as.data.frame()
table_panels_row_pilot = unique(table_panels_row_pilot)
dim(table_panels_row_pilot)
# 1395  6
table_panels_row_pilot$gelID = as.character(table_panels_row_pilot$gelID)

# Let's take only those that have any of the panels specified in part2 within list_panels
table_a_pilot_ultra = table_a_pilot_ultra %>%
  group_by(gelID) %>%
  mutate(any_panel_in_listpanels_part2 = any_exist(list_panels_part2, panel_list)) %>%
  ungroup() %>%
  as.data.frame()
dim(table_a_pilot_ultra)
# 133 14

# just take the ones including any of the panels suggested
table_a_pilot_ultra = table_a_pilot_ultra %>%
  filter(any_panel_in_listpanels_part2)
dim(table_a_pilot_ultra)
# 14  14

table_a_pilot_ultra = table_a_pilot_ultra[,-14]

table_a_pilot = rbind(table_a_pilot,
                      table_a_pilot_HA,
                      table_a_pilot_ultra)
table_a_pilot = unique(table_a_pilot)
dim(table_a_pilot)
# 425 13

# How many PIDs are in the Pilot?
length(unique(table_a_pilot$plateKey))
# 413
length(unique(table_a_pilot$gelID))
# 413

# Let's select the interesting columns for Table 2
table_a = table_a %>% select(participant_id, plate_key.x, rare_diseases_family_id, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, panel_list)
table_a_pilot = table_a_pilot %>% select(gelID, plateKey, gelFamilyId.x, sex, yearOfBirth, specificDisease, panel_list)
colnames(table_a_pilot) = colnames(table_a)

# Let's merge MAIN and PILOT
panel_a = rbind(table_a,
                table_a_pilot)
panel_a = unique(panel_a)
dim(panel_a)
# 3711  7

length(unique(panel_a$participant_id))
# 3692
panel_a$panel = rep("A", length(panel_a$participant_id))

# Enrich it with family-history table
panel_a = left_join(panel_a,
                    is_fami %>% select(rare_diseases_family_id, is_fam),
                    by = "rare_diseases_family_id")
panel_a = left_join(panel_a,
                    is_fami_pilot %>% select(gelFamilyId.x, is_fam_pilot),
                    by = c("rare_diseases_family_id" = "gelFamilyId.x"))

panel_a %>% filter(is_fam) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 1196
panel_a %>% filter(is_fam_pilot) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 140

# per normalised spec disease
panel_a %>% filter(is_fam, normalised_specific_disease %in% "Hereditary ataxia") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 352
panel_a %>% filter(is_fam_pilot, normalised_specific_disease %in% "Hereditary ataxia") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 51

panel_a %>% filter(is_fam, normalised_specific_disease %in% "Hereditary spastic paraplegia") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 182
panel_a %>% filter(is_fam_pilot, normalised_specific_disease %in% "Hereditary spastic paraplegia") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 39

panel_a %>% filter(is_fam, normalised_specific_disease %in% "Early onset and familial Parkinson's Disease") %>% select(participant_id) %>% unique() %>% pull() %>% length()
#
panel_a %>% filter(is_fam_pilot, normalised_specific_disease %in% "Early onset and familial Parkinson's Disease") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 5

panel_a %>% filter(is_fam, grepl("Complex Parkinsonism", normalised_specific_disease)) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 28
panel_a %>% filter(is_fam_pilot, grepl("Complex Parkinsonism", normalised_specific_disease)) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 3

panel_a %>% filter(is_fam, normalised_specific_disease %in% "Early onset dystonia") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 89
panel_a %>% filter(is_fam_pilot, normalised_specific_disease %in% "Early onset dystonia") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 15

panel_a %>% filter(is_fam, normalised_specific_disease %in% "Early onset dementia") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 88
panel_a %>% filter(is_fam_pilot, normalised_specific_disease %in% "Early onset dementia") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 0

panel_a %>% filter(is_fam, normalised_specific_disease %in% "Amyotrophic lateral sclerosis or motor neuron disease") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 19
panel_a %>% filter(is_fam_pilot, normalised_specific_disease %in% "Amyotrophic lateral sclerosis or motor neuron disease") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 0

panel_a %>% filter(is_fam, normalised_specific_disease %in% "Charcot-Marie-Tooth disease") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 242
panel_a %>% filter(is_fam_pilot, normalised_specific_disease %in% "Charcot-Marie-Tooth disease") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 36

panel_a %>% filter(is_fam, normalised_specific_disease %in% "Ultra-rare undescribed monogenic disorders") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 19
panel_a %>% filter(is_fam_pilot, normalised_specific_disease %in% "Ultra-rare undescribed monogenic disorders") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 0

################################################################################################################################################################
# TABLE B - Complex ID
# We need to take the list of PIDs from `list_2459_PIDs_ID_and_others_as_panels.txt`

# load list of 2449 PIDs 
l_complex_ID_group2 = read.table("./list_2743_PIDs_ID_and_others_as_panels.txt", stringsAsFactors = F)
l_complex_ID_group2 = l_complex_ID_group2$V1
length(l_complex_ID_group2)
# 2743

table_b = table_diseases %>%
  filter(participant_id %in% l_complex_ID_group2)
dim(table_b)
# 3132  21

length(unique(table_b$plate_key.x))
# 2743

length(unique(table_b$participant_id))
# 2743

panel_b = table_b %>% select(participant_id, plate_key.x, rare_diseases_family_id, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, panel_list)
panel_b$panel = rep("B", length(panel_b$participant_id))
panel_b = unique(panel_b)
dim(panel_b)
# 2851  8

length(unique(panel_b$participant_id))
# 2743


panel_b = left_join(panel_b,
                    is_fami %>% select(rare_diseases_family_id,is_fam),
                    by = "rare_diseases_family_id")
panel_b = left_join(panel_b,
                    is_fami_pilot %>% select(gelFamilyId.x,is_fam_pilot),
                    by = c("rare_diseases_family_id" = "gelFamilyId.x"))

panel_b %>% filter(is_fam) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 528
panel_b %>% filter(is_fam_pilot) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 0


# Analysing numbers for Panel B
# Define here the list of panels we want to have within
list_panels = c("Genetic epilepsy syndromes", " Genetic epilepsy syndromes",
                "Congenital muscular dystrophy", " Congenital muscular dystrophy",
                "Hereditary ataxia"," Hereditary ataxia",
                " Hereditary spastic paraplegia","Hereditary spastic paraplegia",
                " Mitochondrial disorders","Mitochondrial disorders",
                " Inherited white matter disorders","Inherited white matter disorders",
                "Optic neuropathy", " Optic neuropathy",
                " Brain channelopathy","Brain channelopathy")

# HPO <-> Panels
seizures = c("Genetic epilepsy syndromes", "Epileptic encephalopathy")
dystonia = c("Early onset dystonia")
ataxia = c("Hereditary ataxia", "Brain channelopathy")
spastic_paraplegia = c("Hereditary spastic paraplegia")
optic_neuropathy_or_retinopathy = c("Optic neuropathy")
white_matter_abnormalities = c("Inherited white matter disorders")
muscular_weakness_hypotonia = c("Congenital muscular dystrophy")

# Split panels split by commas
panel_b_panels = panel_b %>% 
  mutate(panels = strsplit(as.character(panel_list), ",")) %>%
  unnest(panels) %>%
  as.data.frame()
panel_b_panels = unique(panel_b_panels)
dim(panel_b_panels)
# 19764  9
panel_b_panels$participant_id = as.character(panel_b_panels$participant_id)
panel_b_panels = panel_b_panels[,-7]
panel_b_panels = unique(panel_b_panels)
dim(panel_b_panels)
# 19764  8

# Let's pull out PIDs having a panel 
l_pids_ID = panel_b_panels %>% filter(grepl("Intellectual disability", panels)) %>% select(participant_id) %>% unique() %>% pull() 

l_pid_dystonia = panel_b_panels %>% filter(grepl(dystonia, panels)) %>% select(participant_id) %>% unique() %>% pull() 

l_pid_seizures_part1 = panel_b_panels %>% filter(grepl("Genetic epilepsy syndromes", panels)) %>% select(participant_id) %>% unique() %>% pull() 
l_pid_seizures_part2 = panel_b_panels %>% filter(grepl("Epileptic encephalopathy", panels)) %>% select(participant_id) %>% unique() %>% pull() 
l_pid_seizures = unique(c(l_pid_seizures_part1,
                          l_pid_seizures_part2))

l_pid_ataxia_part1 = panel_b_panels %>% filter(grepl("Hereditary ataxia", panels)) %>% select(participant_id) %>% unique() %>% pull() 
l_pid_ataxia_part2 = panel_b_panels %>% filter(grepl("Brain channelopathy", panels)) %>% select(participant_id) %>% unique() %>% pull() 
l_pid_ataxia = unique(c(l_pid_ataxia_part1,
                        l_pid_ataxia_part2))

l_pid_spastic_paraplegia = panel_b_panels %>% filter(grepl(spastic_paraplegia, panels)) %>% select(participant_id) %>% unique() %>% pull() 
l_pid_optic_neuro_or_retino = panel_b_panels %>% filter(grepl(optic_neuropathy_or_retinopathy, panels)) %>% select(participant_id) %>% unique() %>% pull() 
l_pid_white = panel_b_panels %>% filter(grepl(white_matter_abnormalities, panels)) %>% select(participant_id) %>% unique() %>% pull() 
l_pid_muscular_hypo = panel_b_panels %>% filter(grepl(muscular_weakness_hypotonia, panels)) %>% select(participant_id) %>% unique() %>% pull() 

# NUmber of PIDs having ID and ONLY seizures
l_pid_seizures_ONLY = setdiff(l_pid_seizures, 
                              c(l_pid_ataxia, l_pid_dystonia, l_pid_spastic_paraplegia, l_pid_optic_neuro_or_retino, l_pid_white, l_pid_muscular_hypo))
length(unique(l_pid_seizures_ONLY))
# 1048

length(unique(intersect(l_pid_seizures_ONLY, l_pids_ID)))
# 1048

l_pid_dystonia_ONLY = setdiff(l_pid_dystonia, 
                              c(l_pid_ataxia, l_pid_seizures, l_pid_spastic_paraplegia, l_pid_optic_neuro_or_retino, l_pid_white, l_pid_muscular_hypo))
length(unique(l_pid_dystonia_ONLY))
# 22
length(unique(intersect(l_pid_dystonia_ONLY, l_pids_ID)))
# 22

l_pid_ataxia_ONLY = setdiff(l_pid_ataxia, 
                            c(l_pid_dystonia, l_pid_seizures, l_pid_spastic_paraplegia, l_pid_optic_neuro_or_retino, l_pid_white, l_pid_muscular_hypo))
length(unique(l_pid_ataxia_ONLY))
# 116

length(unique(intersect(l_pid_ataxia_ONLY, l_pids_ID)))
# 116

l_pid_spastic_paraplegia_ONLY = setdiff(l_pid_spastic_paraplegia, 
                                        c(l_pid_dystonia, l_pid_seizures, l_pid_ataxia, l_pid_optic_neuro_or_retino, l_pid_white, l_pid_muscular_hypo))
length(unique(l_pid_spastic_paraplegia_ONLY))
# 76
length(unique(intersect(l_pid_spastic_paraplegia_ONLY, l_pids_ID)))
# 76

l_pid_optic_neuro_or_retino_ONLY = setdiff(l_pid_optic_neuro_or_retino, 
                                           c(l_pid_dystonia, l_pid_seizures, l_pid_ataxia, l_pid_spastic_paraplegia, l_pid_white, l_pid_muscular_hypo))
length(unique(l_pid_optic_neuro_or_retino_ONLY))
# 10
length(unique(intersect(l_pid_optic_neuro_or_retino_ONLY, l_pids_ID)))
# 10

l_pid_white_ONLY = setdiff(l_pid_white, 
                           c(l_pid_dystonia, l_pid_seizures, l_pid_ataxia, l_pid_spastic_paraplegia, l_pid_optic_neuro_or_retino, l_pid_muscular_hypo))
length(unique(l_pid_white_ONLY))
# 91
length(unique(intersect(l_pid_white_ONLY, l_pids_ID)))
# 91

l_pid_muscular_hypo_ONLY = setdiff(l_pid_muscular_hypo, 
                                   c(l_pid_dystonia, l_pid_seizures, l_pid_ataxia, l_pid_spastic_paraplegia, l_pid_optic_neuro_or_retino, l_pid_white))
length(unique(l_pid_muscular_hypo_ONLY))
# 117
length(unique(intersect(l_pid_muscular_hypo_ONLY, l_pids_ID)))
# 117

# ID and 2 of the above
#l_pid_seizures_ONLY
#l_pid_dystonia_ONLY
#l_pid_ataxia_ONLY
#l_pid_spastic_paraplegia_ONLY
#l_pid_optic_neuro_or_retino_ONLY
#l_pid_white_ONLY
#l_pid_muscular_hypo_ONLY
coctail_panels = c(" Genetic epilepsy syndromes", " Epileptic encephalopathy",
                   " Early onset dystonia",
                   " Hereditary ataxia", " Brain channelopathy",
                   " Hereditary spastic paraplegia",
                   " Optic neuropathy",
                   " Inherited white matter disorders",
                   " Congenital muscular dystrophy")

# example: FID = 111003831 which has ID + HSP + Inherited white mmater disorders + 
panel_b = panel_b %>%
  group_by(participant_id) %>%
  mutate(ID_and_one = ifelse(length(intersect(coctail_panels, unlist(strsplit(panel_list, ",")))) >= 1, TRUE, FALSE)) %>%
  mutate(ID_and_two = ifelse(length(intersect(coctail_panels, unlist(strsplit(panel_list, ",")))) >= 2, TRUE, FALSE)) %>%
  mutate(ID_and_three = ifelse(length(intersect(coctail_panels, unlist(strsplit(panel_list, ",")))) >= 3, TRUE, FALSE)) %>%
  mutate(ID_and_four = ifelse(length(intersect(coctail_panels, unlist(strsplit(panel_list, ",")))) >= 4, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()

length(which(panel_b$ID_and_one))
# 2654
length(which(panel_b$ID_and_two))
# 1022
length(which(panel_b$ID_and_three))
# 615
length(which(panel_b$ID_and_four))
# 317

################################################################################################################################################################
# TABLE C
# patients presenting with intellectual disability and or a neuromuscular phenotype were analysed for DMPK

# MAIN
# Arpil 2021: removing `intellectual disability`
table_c = table_diseases %>%
  filter(normalised_specific_disease %in% c("Congenital muscular dystrophy",
                                            "Congenital myopathy",
                                            "Skeletal Muscle Channelopathies",
                                            "Distal myopathies"))
dim(table_c)
# 805  21

# PILOT
table_c_pilot = table_diseases_pilot %>%
  filter(specificDisease %in% c("Congenital muscular dystrophy",
                                "Congenital myopathy",
                                "Skeletal Muscle Channelopathies",
                                "Distal myopathies"))
dim(table_c_pilot)
# 81  13

table_c = table_c %>% select(participant_id, plate_key.x, rare_diseases_family_id, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, panel_list)
table_c_pilot = table_c_pilot %>% select(gelID, plateKey, gelFamilyId.x, sex, yearOfBirth, specificDisease, panel_list)
colnames(table_c_pilot) = colnames(table_c)

panel_c = rbind(table_c,
                table_c_pilot)
panel_c = unique(panel_c)
dim(panel_c)
# 861  7

panel_c$panel = rep("C", length(panel_c$participant_id))

# How many PIDs
length(unique(panel_c$participant_id))
# 860
length(unique(panel_c$plate_key.x))
# 860

panel_c = left_join(panel_c,
                    is_fami %>% select(rare_diseases_family_id, is_fam),
                    by = "rare_diseases_family_id")
panel_c = left_join(panel_c,
                    is_fami_pilot %>% select(gelFamilyId.x, is_fam_pilot),
                    by = c("rare_diseases_family_id" = "gelFamilyId.x"))

panel_c %>% filter(is_fam) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 193
panel_c %>% filter(is_fam_pilot) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 27

# Per norm spec disease
# Congenital myopathy
panel_c %>% filter(is_fam, normalised_specific_disease %in% "Congenital myopathy") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 101
panel_c %>% filter(is_fam_pilot, normalised_specific_disease %in% "Congenital myopathy") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 15

#Distal myopathies
panel_c %>% filter(is_fam, normalised_specific_disease %in% "Distal myopathies") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 46
panel_c %>% filter(is_fam_pilot, normalised_specific_disease %in% "Distal myopathies") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 6

#Congenital muscular dystrophy
panel_c %>% filter(is_fam, normalised_specific_disease %in% "Congenital muscular dystrophy") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 20
panel_c %>% filter(is_fam_pilot, normalised_specific_disease %in% "Congenital muscular dystrophy") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 4

#Skeletal muscle channelopathy
panel_c %>% filter(is_fam, normalised_specific_disease %in% "Skeletal Muscle Channelopathies") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 26
panel_c %>% filter(is_fam_pilot, normalised_specific_disease %in% "Skeletal Muscle Channelopathies") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 3


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

table_d = table_d %>% select(participant_id, plate_key.x, rare_diseases_family_id, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, panel_list)
table_d_pilot = table_d_pilot %>% select(gelID, plateKey, gelFamilyId.x, sex, yearOfBirth, specificDisease, panel_list)
colnames(table_d_pilot) = colnames(table_d)

panel_d = rbind(table_d,
                table_d_pilot)
panel_d = unique(panel_d)
dim(panel_d)
# 6731  7

panel_d$panel = rep("D", length(panel_d$participant_id))

# PIDs?
length(unique(panel_d$participant_id))
# 6731


panel_d = left_join(panel_d,
                    is_fami %>% select(rare_diseases_family_id,is_fam),
                    by = "rare_diseases_family_id")
panel_d = left_join(panel_d,
                    is_fami_pilot %>% select(gelFamilyId.x,is_fam_pilot),
                    by = c("rare_diseases_family_id" = "gelFamilyId.x"))

panel_d %>% filter(is_fam) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 1481
panel_d %>% filter(is_fam_pilot) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 55

# Let's merge all panels in the same table
panel_merged = rbind(panel_a,
                     panel_b,
                     panel_c,
                     panel_d)
panel_merged = unique(panel_merged)
dim(panel_merged)
# 14154  10

panel_merged %>% filter(is_fam) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 2921 
panel_merged %>% filter(is_fam_pilot) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 218

panel_merged = panel_merged %>% group_by(participant_id) %>% mutate(age = 2020 - year_of_birth) %>% ungroup() %>% as.data.frame()

# Age and sex
summary(panel_merged$age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.00    9.00   16.00   25.69   39.00  101.00 
panel_merged %>% select(participant_id, participant_phenotypic_sex) %>% unique() %>% select(participant_phenotypic_sex) %>% table()
#female Female   male   Male 
#302   4652    348   6329

panel_merged = rbind(panel_a,
                     panel_c,
                     panel_d)
panel_merged = unique(panel_merged)
dim(panel_merged)
# 18034  8


# PIDs?
length(unique(panel_merged$participant_id))
# 11631 with A,B,C,D
# 11266 with A,C,D

# Calculate for each disease in Table2 summary of age and gender distribution
l_diseases_table2 = c("Amyotrophic lateral sclerosis or motor neuron disease", 
                      "Charcot-Marie-Tooth disease",
                      "Early onset dementia", 
                      "Early onset dystonia", 
                      "Complex Parkinsonism", 
                      "Hereditary ataxia", 
                      "Hereditary spastic paraplegia",
                      "'Early onset and familial Parkinson''s Disease'",
                      "Intellectual disability",
                      "Kabuki syndrome",
                      "Congenital muscular dystrophy",
                      "Congenital myopathy",
                      "Skeletal Muscle Channelopathies",
                      "Distal myopathies",
                      "Ultra-rare undescribed monogenic disorders",
                      "Unknown disorder", 
                      "All recognised syndromes and those with suggestive features")

panel_merged = panel_merged %>%
  group_by(participant_id) %>%
  mutate(age = 2020 - year_of_birth) %>%
  ungroup() %>%
  as.data.frame()

for(i in 1:length(l_diseases_table2)){
  print(l_diseases_table2[i])
  
  l_pid_disease = panel_merged %>% filter(normalised_specific_disease %in% l_diseases_table2[i]) %>% select(participant_id) %>% unique() %>% pull()
  print(length(l_pid_disease))
  
  panel_merged %>% filter(participant_id %in% l_pid_disease) %>% select(age) %>% summary() %>% print()
  
  # Panel A
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "A") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "B") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "C") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "D") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
  
  # families
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "A") %>% select(rare_diseases_family_id) %>% unique() %>% pull() %>% length() %>% print()
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "B") %>% select(rare_diseases_family_id) %>% unique() %>% pull() %>% length() %>% print()
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "C") %>% select(rare_diseases_family_id) %>% unique() %>% pull() %>% length() %>% print()
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "D") %>% select(rare_diseases_family_id) %>% unique() %>% pull() %>% length() %>% print()
  
  
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "A") %>% select(participant_phenotypic_sex, participant_id) %>% unique() %>% select(participant_phenotypic_sex) %>% table() %>% print()
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "B") %>% select(participant_phenotypic_sex, participant_id) %>% unique() %>% select(participant_phenotypic_sex) %>% table() %>% print()
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "C") %>% select(participant_phenotypic_sex, participant_id) %>% unique() %>% select(participant_phenotypic_sex) %>% table() %>% print()
  panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "D") %>% select(participant_phenotypic_sex, participant_id) %>% unique() %>% select(participant_phenotypic_sex) %>% table() %>% print()
}
#complex parkin
l_complex = panel_merged %>% filter(grepl("[Cc]omplex [Pp]arkin",normalised_specific_disease, ignore.case = T)) %>% select(participant_id) %>% unique() %>% pull()
print(length(l_complex))

panel_merged %>% filter(grepl("[Cc]omplex [Pp]arkin",normalised_specific_disease, ignore.case = T)) %>% select(age) %>% summary() %>% print()

panel_merged %>% filter(grepl("[Cc]omplex [Pp]arkin",normalised_specific_disease, ignore.case = T)) %>% select(participant_phenotypic_sex, participant_id) %>% unique() %>% select(participant_phenotypic_sex) %>% table() %>% print()

panel_merged %>% filter(grepl("[Cc]omplex [Pp]arkin",normalised_specific_disease, ignore.case = T), panel %in% "A") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
panel_merged %>% filter(grepl("[Cc]omplex [Pp]arkin",normalised_specific_disease, ignore.case = T), panel %in% "B") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
panel_merged %>% filter(grepl("[Cc]omplex [Pp]arkin",normalised_specific_disease, ignore.case = T), panel %in% "C") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
panel_merged %>% filter(grepl("[Cc]omplex [Pp]arkin",normalised_specific_disease, ignore.case = T), panel %in% "D") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()


# Ultra-rare in MAIN
l_pid_disease = panel_merged %>% filter(normalised_specific_disease %in% "Ultra-rare undescribed monogenic disorders") %>% select(participant_id) %>% unique() %>% pull()
print(length(l_pid_disease))

panel_merged %>% filter(participant_id %in% l_pid_disease) %>% select(age) %>% summary() %>% print()

panel_merged %>% filter(participant_id %in% l_pid_disease) %>% select(participant_phenotypic_sex, participant_id) %>% unique() %>% select(participant_phenotypic_sex) %>% table() %>% print()

# Ultra-rare in PILOT
#"Unknown disorder",
#"All recognised syndromes and those with suggestive features"
l_pid_disease = panel_merged %>% 
  filter(normalised_specific_disease %in% c("Unknown disorder", "All recognised syndromes and those with suggestive features")) %>% 
  select(participant_id) %>% unique() %>% pull()
print(length(l_pid_disease))

panel_merged %>% filter(participant_id %in% l_pid_disease) %>% select(age) %>% summary() %>% print()

panel_merged %>% filter(participant_id %in% l_pid_disease) %>% select(participant_phenotypic_sex, participant_id) %>% unique() %>% select(participant_phenotypic_sex) %>% table() %>% print()

panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "A") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "B") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "C") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()
panel_merged %>% filter(participant_id %in% l_pid_disease, panel %in% "D") %>% select(participant_id) %>% unique() %>% pull() %>% length() %>% print()


panel_merged_without_B = panel_merged %>%
  filter(!panel %in% "B")

length(unique(panel_merged_without_B$participant_id))
# 11266

# Age of onset
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/raw/rare_diseases_participant_dise_2020-12-30_12-10-36.tsv",
                     sep = "\t",
                     stringsAsFactors = F)
dim(clin_data)
# 39676  11

# Georgia Chan worked across HES tables to fish missing data re age on onset
georgia = read.csv("~/Documents/STRs/PAPERS/VALIDATION_PAPER/LANCET/APPEAL/Neuro_DON_georgiaChan.csv",
                   stringsAsFactors = F,
                   header = T)
dim(georgia)
# 1228 2

# We need to know the age of the PIDs from Georgia, to define the age of onset (from the date of onset)
main_clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V12_and_Pilot_programmes.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(main_clin_data)
# 2472865  26

pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/data_freeze_Pilot_LK_RESEARCH/pedigree.csv",
                           stringsAsFactors = F,
                           header = T)
dim(pilot_clin_data)
# 17258  34

georgia = left_join(georgia,
                    main_clin_data %>% select(participant_id, year_of_birth),
                    by = "participant_id")

georgia = georgia %>%
  group_by(participant_id) %>%
  mutate(year_onset  = strsplit(date_of_onset, '-')[[1]][1]) %>%
  ungroup() %>%
  as.data.frame() %>%
  unique()

georgia = georgia %>%
  group_by(participant_id) %>%
  mutate(age_of_onset  = as.integer(year_onset) - as.integer(year_of_birth)) %>%
  ungroup() %>%
  as.data.frame() %>%
  unique()

# also include in georgia the normalised disease
georgia = left_join(georgia,
                    clin_data %>% select(participant_id, normalised_specific_disease),
                    by = "participant_id")
  

l_pid_all_panels = unique(c(panel_a$participant_id,
                            panel_b$participant_id,
                            panel_c$participant_id,
                            panel_d$participant_id))
clin_data = clin_data %>% filter(participant_id %in% l_pid_all_panels) %>% select(participant_id, normalised_age_of_onset, normalised_specific_disease) %>% unique()

# Overall
mean(c(clin_data$normalised_age_of_onset,georgia$age_of_onset), na.rm = T)
# 12.42
sd(c(clin_data$normalised_age_of_onset, georgia$age_of_onset), na.rm = T)
# 20.22

# Age of onset Panel A,B,C,D
l_panel_a_age = clin_data %>% filter(participant_id %in% panel_a$participant_id, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull()
mean(c(l_panel_a_age, georgia %>% filter(participant_id %in% panel_a$participant_id) %>% select(age_of_onset) %>% pull()), na.rm = T)
# 34.84
sd(c(l_panel_a_age, georgia %>% filter(participant_id %in% panel_a$participant_id) %>% select(age_of_onset) %>% pull()), na.rm = T)
# 21.39

l_panel_b_age = clin_data %>% filter(participant_id %in% panel_b$participant_id, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() 
mean(c(l_panel_b_age, georgia %>% filter(participant_id %in% panel_b$participant_id) %>% select(age_of_onset) %>% pull()), na.rm = T)
# 1.60 
sd(c(l_panel_b_age, georgia %>% filter(participant_id %in% panel_b$participant_id) %>% select(age_of_onset) %>% pull()), na.rm = T)
# 5.28

l_panel_c_age = clin_data %>% filter(participant_id %in% panel_c$participant_id, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() 
mean(c(l_panel_c_age, georgia %>% filter(participant_id %in% panel_c$participant_id) %>% select(age_of_onset) %>% pull()), na.rm = T)
# 17.75
sd(c(l_panel_c_age, georgia %>% filter(participant_id %in% panel_c$participant_id) %>% select(age_of_onset) %>% pull()), na.rm = T)
# 21.06

l_panel_d_age = clin_data %>% filter(participant_id %in% panel_d$participant_id, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull()
mean(c(l_panel_d_age, georgia %>% filter(participant_id %in% panel_d$participant_id) %>% select(age_of_onset) %>% pull()), na.rm = T)
# 1.06
sd(c(l_panel_d_age, georgia %>% filter(participant_id %in% panel_d$participant_id) %>% select(age_of_onset) %>% pull()), na.rm = T)
# 3.12

# Age of onset on confirmed platekeys across panels
confirmed_panelA = c("112008033","112005899","112005900","115008391","115014127","210012854","115008056","115000968","115006565","116002069","115004635","113004320","115000536","115017486","115017493","115011422","115016317","115006729","113002006","119001845","115002538","115011586","115011457","115007340","115004045","122001169","115004196","210013915","50003892","115012847","115011956","115001026","113001410","50001165","50002396","112003594","113003032","113002194","120000698","111001710","111001708","124000946","122000175","115005431","118002342","118001261")
clin_data %>% filter(participant_id %in% confirmed_panelA, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% confirmed_panelA, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

# Panel A and Hereditary ataxia
clin_data %>% filter(participant_id %in% c("115000968","116002069","115004635","113004320","115017486","115016317","113002006","119001845","115002538","50003892","115012847","115011956","115001026","50001165","113003032","113002194","120000698","124000946","122000175"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("115000968","116002069","115004635","113004320","115017486","115016317","113002006","119001845","115002538","50003892","115012847","115011956","115001026","50001165","113003032","113002194","120000698","124000946","122000175"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd() 

# Panel A and Hereditary spastic paraplegia
clin_data %>% filter(participant_id %in% c("113004320", "50003892", "50002396"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("113004320", "50003892", "50002396"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

georgia %>% filter(participant_id %in% c("113004320", "50003892", "50002396")) %>% select(age_of_onset) %>% pull() %>% mean() 
georgia %>% filter(participant_id %in% c("113004320", "50003892", "50002396")) %>% select(age_of_onset) %>% pull() %>% sd()


# Panel A and Early onset and familial Parkinson's Disease
clin_data %>% filter(participant_id %in% c("115000536","115004196"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("115000536","115004196"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

georgia %>% filter(participant_id %in% c("115000536","115004196")) %>% select(age_of_onset) %>% pull() %>% mean() 
georgia %>% filter(participant_id %in% c("115000536","115004196")) %>% select(age_of_onset) %>% pull() %>% sd()


# Panel A and Complex Parkinsonism (includes pallido-pyramidal syndromes)
clin_data %>% filter(participant_id %in% c("115006729","112003594"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("115006729","112003594"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

# Panel A and Early onset dementia
clin_data %>% filter(participant_id %in% c("115005431","115011457","122001169","210013915"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("115005431","115011457","122001169","210013915"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

  # Panel A and Amyotrophic lateral sclerosis or motor neuron disease
clin_data %>% filter(participant_id %in% c("112005899","112005900","115008391","115006565","115011422","115011586","115007340","115004045"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("112005899","112005900","115008391","115006565","115011422","115011586","115007340","115004045"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

# Panel A and Charcot-Marie-Tooth disease
clin_data %>% filter(participant_id %in% c("112008033","115014127","210012854","115008056"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("112008033","115014127","210012854","115008056"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

# Panel A and Ultra-rare undescribed monogenic disorders
clin_data %>% filter(participant_id %in% c("113001410","111001710","111001708"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("113001410","111001710","111001708"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

confirmed_panelB = c("111001410","116000121","115001700","118002342","115005431","210013360","115014499")
clin_data %>% filter(participant_id %in% confirmed_panelB, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% confirmed_panelB, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

confirmed_panelC = c("122007152","122007274","210017355","210013126","210013125")
clin_data %>% filter(participant_id %in% confirmed_panelC, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% confirmed_panelC, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

# Panel C and Congenital myopathy
clin_data %>% filter(participant_id %in% "210017355", !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% "210017355", !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

# Panel C and Distal myopathies
clin_data %>% filter(participant_id %in% c("122007152", "122007274"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("122007152", "122007274"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd() 


# Panel C and Congenital muscular dystrophy
clin_data %>% filter(participant_id %in% c("210013126","210013125"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% c("210013126","210013125"), !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

georgia %>% filter(participant_id %in% "210013126") %>% select(age_of_onset) %>% pull() %>% mean() 
georgia %>% filter(participant_id %in% "210013126") %>% select(age_of_onset) %>% pull() %>% sd()

# Panel C and Skeletal muscle channelopathy

confirmed_panelD = c("118002794","117000919","112002287","122005899","112001315","112001329","112001252","115005821","116001580","116000367")
clin_data %>% filter(participant_id %in% confirmed_panelD, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% confirmed_panelD, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()

all_confirmed_panels = unique(c(confirmed_panelA,
                         confirmed_panelB,
                         confirmed_panelC,
                         confirmed_panelD))
length(all_confirmed_panels)
# 66

clin_data %>% filter(participant_id %in% all_confirmed_panels, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() 
clin_data %>% filter(participant_id %in% all_confirmed_panels, !is.na(normalised_age_of_onset)) %>% select(normalised_age_of_onset) %>% pull() %>% sd()


# Age of onset across different diseases in diff panels
for (i in 1:length(l_diseases_table2)){
  print(l_diseases_table2[i])
  print("Panel A")
  aux_a = clin_data %>% filter(participant_id %in% panel_a$participant_id, !is.na(normalised_age_of_onset), normalised_specific_disease %in% l_diseases_table2[i]) %>% select(normalised_age_of_onset) %>% pull()
  print(mean(c(aux_a, 
               georgia %>% filter(participant_id %in% panel_a$participant_id, normalised_specific_disease %in% l_diseases_table2[i]) %>% select(age_of_onset) %>% pull()), na.rm = T))
  print(sd(c(aux_a, 
               georgia %>% filter(participant_id %in% panel_a$participant_id, normalised_specific_disease %in% l_diseases_table2[i]) %>% select(age_of_onset) %>% pull()), na.rm = T))
  
  print("Panel C")
  aux_b = clin_data %>% filter(participant_id %in% panel_c$participant_id, !is.na(normalised_age_of_onset), normalised_specific_disease %in% l_diseases_table2[i]) %>% select(normalised_age_of_onset) %>% pull()
  print(mean(c(aux_b, 
               georgia %>% filter(participant_id %in% panel_c$participant_id, normalised_specific_disease %in% l_diseases_table2[i]) %>% select(age_of_onset) %>% pull()), na.rm = T))
  print(sd(c(aux_b, 
               georgia %>% filter(participant_id %in% panel_c$participant_id, normalised_specific_disease %in% l_diseases_table2[i]) %>% select(age_of_onset) %>% pull()), na.rm = T))
  
  print("Panel D")
  aux_d = clin_data %>% filter(participant_id %in% panel_d$participant_id, !is.na(normalised_age_of_onset), normalised_specific_disease %in% l_diseases_table2[i]) %>% select(normalised_age_of_onset) %>% pull()
  print(mean(c(aux_d,
             georgia %>% filter(participant_id %in% panel_d$participant_id, normalised_specific_disease %in% l_diseases_table2[i]) %>% select(age_of_onset) %>% pull()), na.rm = T))
  print(sd(c(aux_d,
               georgia %>% filter(participant_id %in% panel_d$participant_id, normalised_specific_disease %in% l_diseases_table2[i]) %>% select(age_of_onset) %>% pull()), na.rm = T))
  
}

print("Panel A")
clin_data %>% filter(participant_id %in% panel_a$participant_id, !is.na(normalised_age_of_onset), grepl("Complex parkin", normalised_specific_disease, ignore.case = T)) %>% select(normalised_age_of_onset) %>% pull() %>% mean() %>% print()
clin_data %>% filter(participant_id %in% panel_a$participant_id, !is.na(normalised_age_of_onset), grepl("Complex parkin", normalised_specific_disease, ignore.case = T)) %>% select(normalised_age_of_onset) %>% pull() %>% sd() %>% print()

##############

to_check_ultra_platekeys = table_diseases %>% filter(normalised_specific_disease %in% "Ultra-rare undescribed monogenic disorders") %>% select(plate_key.x) %>% unique() %>% pull() 

# Check 1517 ultra-rare have expansions ?
# We need to load here the repeat-size estimations for all them - EHv255
repeats_table_main = read.csv("~/Documents/STRs/data/research/EH_2.5.5_research_August2019/EH_output_v2.5.5_research_August_2019/merged_genotypeUpdated/merged_loci_86457_research_genomes_new_loci_EHv2.5.5_summer2019_removingListVcfFiles.tsv",
                              stringsAsFactors = F,
                              header = T,
                              sep = "\t")
dim(repeats_table_main)
# 3983  11

# recode repeats_table_main: FMR1 for b37 is `FMR1` and `FMR1_CGG` for b38, let's recode
repeats_table_main$gene = recode(repeats_table_main$gene,
                                 "FMR1" = "FMR1_CGG")
# Load the pathogenic threshold for the loci
gene_pathogenic_threshold = read.csv("~/git/analysing_STRs/threshold_smallest_pathogenic_reported.txt",
                                     sep = "\t",
                                     stringsAsFactors = F)

# Now, we want to see how many of them have an expansion on any of the genes in `l_genes_tableA`
l_genes_tableA = c("AR_CAG", "ATN1_CAG", "ATXN1_CAG", "ATXN2_CAG", "ATXN3_CAG", "ATXN7_CAG", "CACNA1A_CAG", "C9orf72_GGGGCC", "FXN_GAA", "HTT_CAG", "TBP_CAG", "FMR1_CGG")

expanded_table_main = data.frame()
for (i in 1:length(l_genes_tableA)){
  locus_name = l_genes_tableA[i]
  patho_cutoff = gene_pathogenic_threshold %>% 
    filter(locus %in% locus_name) %>%
    select(threshold) %>%
    pull()
  
  print(locus_name)
  print(patho_cutoff)
  
  expanded_table_main = rbind(expanded_table_main,
                              repeats_table_main %>% 
                                filter(gene %in% locus_name, allele >= patho_cutoff) %>%
                                select(gene, allele, Repeat_Motif, num_samples, list_samples))
  
}
dim(expanded_table_main)
# 423  5

# Now, we want to see how many of them have an expansion on any of the genes in `l_genes_tableA` - but for Pilot data
expanded_table_pilot = data.frame()
for (i in 1:length(l_genes_tableA)){
  locus_name = l_genes_tableA[i]
  patho_cutoff = gene_pathogenic_threshold %>% 
    filter(locus %in% locus_name) %>%
    select(threshold) %>%
    pull()
  
  print(locus_name)
  print(patho_cutoff)
  
  expanded_table_pilot = rbind(expanded_table_pilot,
                               repeats_table_pilot %>% 
                                 filter(gene %in% locus_name, allele >= patho_cutoff) %>%
                                 select(gene, allele, Repeat_Motif, num_samples, list_samples))
  
}
dim(expanded_table_pilot)
# 69  5

expanded_table_main_per_locus = data.frame()
index_kutre = 1
for (i in 1:length(expanded_table_main$gene)){
  list_affected_vcf = strsplit(expanded_table_main$list_samples[i], ';')[[1]]
  for (j in 1:length(list_affected_vcf)){
    expanded_table_main_per_locus = rbind(expanded_table_main_per_locus,
                                          expanded_table_main[i,])
    expanded_table_main_per_locus$list_samples[index_kutre] = sub(".vcf", "", sub("EH_", "", list_affected_vcf[j]))
    index_kutre = index_kutre + 1
  }
}
expanded_table_main_per_locus = unique(expanded_table_main_per_locus)
dim(expanded_table_main_per_locus)
# 3718  5

# The same for PILOT
expanded_table_pilot_per_locus = data.frame()
index_kutre = 1
for (i in 1:length(expanded_table_pilot$gene)){
  list_affected_vcf = strsplit(expanded_table_pilot$list_samples[i], ';')[[1]]
  for (j in 1:length(list_affected_vcf)){
    expanded_table_pilot_per_locus = rbind(expanded_table_pilot_per_locus,
                                           expanded_table_pilot[i,])
    expanded_table_pilot_per_locus$list_samples[index_kutre] = sub(".vcf", "", sub("EH_", "", list_affected_vcf[j]))
    index_kutre = index_kutre + 1
  }
}
expanded_table_pilot_per_locus = unique(expanded_table_pilot_per_locus)
dim(expanded_table_pilot_per_locus)
# 176  5

# From the expanded table, let's see how many are in l_platekeys_tableA
expanded_table_main_in_tableA = expanded_table_main_per_locus %>%
  filter(list_samples %in% to_check_ultra_platekeys)
dim(expanded_table_main_in_tableA)
# 65  5

# The same por PILOT
expanded_table_pilot_in_tableA = expanded_table_pilot_per_locus %>%
  filter(list_samples %in% to_check_ultra_platekeys)
dim(expanded_table_pilot_in_tableA)
# 0  5

write.table(expanded_table_main_in_tableA, "./ultra-rare_expanded_MAIN.tsv", quote = F, col.names = T, row.names = F, sep = "\t")
