# Objective: select the platekeys/participant ids corresponding to the list of diseases we consider for the paper
# Both: main and pilot programmes
# libraries
library(dplyr)
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1

setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

# `clinical_data_research_cohort_86457_genomes_withPanels_250919` is a table I generated
clinical_data_research_cohort_86457_genomes_withPanels_250919 = read.csv("~/Documents/STRs/clinical_data/clinical_research_cohort/clinical_data_research_cohort_86457_genomes_withPanels_041019.tsv",
                                                                         sep = "\t",
                                                                         stringsAsFactors = F,
                                                                         header = T)
dim(clinical_data_research_cohort_86457_genomes_withPanels_250919)
# 89305  16

clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V11_and_Pilot_programmes.tsv",
                     stringsAsFactors = F,
                     sep = "\t")
dim(clin_data)
# 2468506  26

list_diseases = c("Intellectual disability",
                  "Amyotrophic lateral sclerosis or motor neuron disease", 
                  "Charcot-Marie-Tooth disease", "Congenital muscular dystrophy",
                  "Congenital myopathy", "Early onset dementia", "Early onset dystonia", 
                  "Distal myopathies", "Complex Parkinsonism", "Hereditary ataxia", 
                  "Hereditary spastic paraplegia", "Skeletal Muscle Channelopathies",
                  "'Early onset and familial Parkinson''s Disease'",
                  "Mitochondrial disorders",
                  "Kabuki syndrome",
                  "Ultra-rare undescribed monogenic disorders")

table_diseases = clin_data %>%
  filter(grepl("Intellectual disability", diseases_list, ignore.case = T) | 
           grepl("Amyotrophic lateral", diseases_list, ignore.case = T) |
           grepl("Charcot-Marie-Tooth disease", diseases_list, ignore.case = T) |
           grepl("Congenital myopathy", diseases_list, ignore.case = T) |
           grepl("Early onset dementia", diseases_list, ignore.case = T) |
           grepl("Congenital muscular dystrophy", diseases_list, ignore.case = T) |
           grepl("Early onset dystonia", diseases_list, ignore.case = T) |
           grepl("Distal myopathies", diseases_list, ignore.case = T) |
           grepl("Complex Parkinsonism", diseases_list, ignore.case = T) |
           grepl("Hereditary ataxia", diseases_list, ignore.case = T) |
           grepl("Hereditary spastic paraplegia", diseases_list, ignore.case = T) |
           grepl("Skeletal Muscle Channelopathies", diseases_list, ignore.case = T) |
           grepl("Early onset and familial Parkinson's Disease", diseases_list, ignore.case = T) |
           #grepl("Mitochondrial disorders", diseases_list, ignore.case = T) |
           grepl("Kabuki syndrome", diseases_list, ignore.case = T) |
           grepl("Ultra-rare undescribed monogenic disorders", diseases_list, ignore.case = T) | 
           grepl("Complex Parkin", diseases_list, ignore.case = T))
dim(table_diseases)           
# 1283953  26 

# Check unique number of platekeys and PIDs
length(unique(table_diseases$platekey))
# 14534
length(unique(table_diseases$participant_id))
# 14421

table_diseases_dedup = table_diseases %>% 
  group_by(participant_id) %>%
  mutate(latest_platekey = max(platekey)) %>%
  ungroup() %>%
  as.data.frame()

table_diseases_dedup = table_diseases_dedup[,-2]
table_diseases_dedup = unique(table_diseases_dedup)
dim(table_diseases_dedup)
# 28300  26

write.table(table_diseases_dedup,
            "./table_diseases_for_table2_Main_and_Pilot_14785_PIDs_all_adults_and_paediatrics_withouth_Mitochondrial.tsv",
            quote = F,
            row.names = F, col.names = T, sep = "\t")

table_diseases <-  clinical_data_research_cohort_86457_genomes_withPanels_250919 %>%
  count(normalised_specific_disease)

count_of_diseases <- clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% count(normalised_specific_disease)

table_diseases <- clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% 
  filter(normalised_specific_disease %in% c("Intellectual disability",
                                            "Amyotrophic lateral sclerosis or motor neuron disease", 
                                            "Charcot-Marie-Tooth disease", "Congenital muscular dystrophy",
                                            "Congenital myopathy", "Early onset dementia", "Early onset dystonia", 
                                            "Distal myopathies", "Complex Parkinsonism", "Hereditary ataxia", 
                                            "Hereditary spastic paraplegia", "Skeletal Muscle Channelopathies",
                                            "'Early onset and familial Parkinson''s Disease'",
                                            "Mitochondrial disorders",
                                            "Kabuki syndrome",
                                            "Ultra-rare undescribed monogenic disorders"))
dim(table_diseases)
# 13237  16

table(table_diseases$normalised_specific_disease)

parkinson_to_enrich = clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% 
  filter(grepl("Complex Parkin", normalised_specific_disease))
dim(parkinson_to_enrich)
# 141  16

table(parkinson_to_enrich$normalised_specific_disease)

table_diseases = rbind(table_diseases, parkinson_to_enrich)
dim(table_diseases)
# 13378  16

write.table(table_diseases[ ,2], row.names = FALSE, col.names = FALSE, quote = FALSE, "./list_PIDs_table_diseases.txt")

# Here we need to deduplicate genomes and take the unique PIDs
table_diseases_dedup = table_diseases %>%
  select(participant_id, plate_key.x, genome_build)

table_diseases_dedup = unique(table_diseases_dedup)
dim(table_diseases_dedup)
# 12978  3

# Create new variable named latest_platekey
table_diseases_dedup =  table_diseases_dedup %>%
  group_by(participant_id) %>%
  mutate(latest_platekey = max(plate_key.x)) %>%
  ungroup() %>%
  as.data.frame()
table_diseases_dedup = table_diseases_dedup %>%
  select(participant_id, latest_platekey, genome_build)
table_diseases_dedup = unique(table_diseases_dedup)
dim(table_diseases_dedup)
# 12905  3

list_duplicated_pid = table_diseases_dedup$participant_id[which(duplicated(table_diseases_dedup$participant_id))]
length(list_duplicated_pid)
# 219

# Take GRCh38 platekey for duplicate PIDs
table_diseases_duplicates = table_diseases_dedup %>%
  filter(participant_id %in% list_duplicated_pid)

table_diseases_duplicates = table_diseases_duplicates %>%
  filter(genome_build %in% "GRCh38")
dim(table_diseases_duplicates)
# 219  3

table_diseases_dedup = table_diseases_dedup %>%
  filter(!participant_id %in% list_duplicated_pid)
dim(table_diseases_dedup)
# 12467  3

# merge
table_diseases_dedup = rbind(table_diseases_dedup,
                             table_diseases_duplicates)

dim(table_diseases_dedup)
# 12686  3

l_platekeys = table_diseases_dedup$latest_platekey
length(l_platekeys)
# 12686

table_diseases_enriched = clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% 
  filter(plate_key.x %in% l_platekeys)
dim(table_diseases_enriched)
# 13868  16

write.table(table_diseases_enriched, file = "table_diseases_enriched_including_skeletalMuscleChan_and_Ultra-rare.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Include or merge `table_diseases_enriched` with PILOT data following the same strategy
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_withPanels_280919.tsv",
                           header = T,
                           sep = "\t",
                           stringsAsFactors = F)
dim(pilot_clin_data)
# 4974  11

# Select individuals/genomes having recruited under any of the following diseases, for the pilot table `specificDisease`
#"Intellectual disability",
#"Amyotrophic lateral sclerosis or motor neuron disease", 
#"Charcot-Marie-Tooth disease", "Congenital muscular dystrophy",
#"Congenital myopathy", "Early onset dementia", "Early onset dystonia", 
#"Distal myopathies", "Complex Parkinsonism", "Hereditary ataxia", 
#"Hereditary spastic paraplegia", "Skeletal Muscle Channelopathies",
#"'Early onset and familial Parkinson''s Disease'"))
table_diseases_pilot = pilot_clin_data %>%
  filter(specificDisease %in% c("Intellectual disability",
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
                                "Early onset and familial Parkinson's Disease",
                                "Unknown disorder",
                                "All recognised syndromes and those with suggestive features"))

View(table(table_diseases_pilot$specificDisease))
# It's ok, we have them all
dim(table_diseases_pilot)
# 660  11

length(unique(table_diseases_pilot$plateKey))
# 645
# 810 (arpil'21)

# There are 645 unique genomes that have been recruited under any of these 13 diseases, some of them (15) have 1 or more above specificDisease
# I won't merge main with pilot since we do have different platekeys and column names.

write.table(table_diseases_pilot, file = "table_diseases_enriched_PILOT_13diseases_22April2021.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Enrich this table with popu  - to take best_guess-predicted_ancestry
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      sep = ",",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59464  36

popu_table_pilot = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            header = T)
dim(popu_table_pilot)
# 4821  44

table_diseases_enriched_popu = left_join(table_diseases_enriched,
                                    popu_table %>% select(ID,best_guess_predicted_ancstry, self_reported),
                                    by = c("plate_key.x" = "ID"))

table_diseases_pilot_popu = left_join(table_diseases_pilot,
                                      popu_table_pilot %>% select(ID, bestGUESS_sub_pop, bestGUESS_super_pop),
                                      by = c("plateKey" = "ID"))

# take reported ancestry
rd_genomes_re = read.table("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_300320.tsv",
                       sep = "\t",
                       stringsAsFactors = FALSE, 
                       header = TRUE)
dim(rd_genomes_re)  
# 1124633  31

table_diseases_enriched_popu = left_join(table_diseases_enriched_popu,
                                         rd_genomes_re %>% select(platekey, participant_ethnic_category),
                                         by = c("plate_key.x" = "platekey"))

table_diseases_enriched_popu = unique(table_diseases_enriched_popu)
dim(table_diseases_enriched_popu)
# 13868 19

write.table(table_diseases_enriched_popu, "table_diseases_enriched_popu_includingSkeletalMuscleChan_and_ultra-rare.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(table_diseases_pilot_popu, "table_diseases_enriched_PILOT_13diseases_enriched_popu.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# Distinguish participants/genomes having Intellectual disability as panel
# Group 1: those having ONLY Intellectual disability in `panel_list`
# Group 2: those having Intellectual disability AND either of the following ones: `Intellectual disability`, `Epileptic encephalopathy`, 
# `Genetic epilepsy syndromes`

# let's make life simple: split into rows panel info
table_panels_row = table_diseases_enriched_popu %>% 
  select(plate_key.x, participant_id, normalised_specific_disease, disease_sub_group, disease_group, panel_list) %>%
  mutate(panels = strsplit(as.character(panel_list), ",")) %>%
  unnest(panels) %>%
  as.data.frame()
table_panels_row = unique(table_panels_row)
dim(table_panels_row)
# 58449  7

table_panels_row$participant_id = as.character(table_panels_row$participant_id)

# I have seen that Intellectual disability in panels appears in two ways:
#`Intellectual disability` and ` Intellectual disability`
# Let's recode
table_panels_row$panels = recode(table_panels_row$panels,
                                 " Intellectual disability" = "Intellectual disability")

table_panels_row = unique(table_panels_row)
dim(table_panels_row)
# 55355  7


# Group 1
l_pid_only_one_panel = table_panels_row %>% 
  group_by(participant_id) %>% 
  filter(n()==1) %>%
  select(participant_id) %>%
  unique() %>%
  pull() %>%
  as.character()
length(l_pid_only_one_panel)  
# 2450

# From 2450 PIDs having only A UNIQUE panel assigned, which ones have been assigned ID
l_pid_ID_group1 = table_panels_row %>%
  filter(grepl("Intellectual disability", panel_list) & participant_id %in% l_pid_only_one_panel) %>%
  select(participant_id) %>%
  unique() %>%
  pull() %>%
  as.character()
length(l_pid_ID_group1)
# 2173

write.table(l_pid_ID_group1, "list_2173_PIDs_only_ID_as_panel_assigned.txt", quote = F, row.names = F, col.names = F)


# Group 2
# Define here the list of panels we want to have within
list_panels = c("Genetic epilepsy syndromes", " Genetic epilepsy syndromes",
                "Congenital muscular dystrophy", " Congenital muscular dystrophy",
                "Hereditary ataxia"," Hereditary ataxia",
                " Hereditary spastic paraplegia","Hereditary spastic paraplegia",
                " Mitochondrial disorders","Mitochondrial disorders",
                " Inherited white matter disorders","Inherited white matter disorders",
                "Optic neuropathy", " Optic neuropathy",
                " Brain channelopathy","Brain channelopathy")




# Function that checks if any of the items in list of characters 1 does exist in list of characters 2
any_exist <- function(list1, list2) {
  for (i in list1){
    if (i %in% list2){
      return(TRUE)
    }
  }
  return(FALSE)
}


l_ID_group2 = table_panels_row %>%
  group_by(participant_id) %>%
  filter(grepl("Intellectual disability", panel_list) & any_exist(list_panels,panels)) %>%
  select(participant_id) %>%
  unique() %>%
  pull() %>%
  as.character()
length(l_ID_group2)
# 2743

# What about PILOT?
table_panels_pilot = table_diseases_pilot %>% 
  select(plateKey, gelID, specificDisease, panel_list) %>%
  mutate(panels = strsplit(as.character(panel_list), ",")) %>%
  unnest(panels) %>%
  as.data.frame()
table_panels_pilot = unique(table_panels_pilot)
dim(table_panels_pilot)
# 1129  5

l_ID_group2_pilot = table_panels_pilot %>%
  group_by(gelID) %>%
  filter(grepl("Intellectual disability", panels) & any_exist(list_panels,panels)) %>%
  select(gelID) %>%
  unique() %>%
  pull() %>%
  as.character()
length(l_ID_group2_pilot)
# 12


write.table(l_ID_group2, "./list_2743_PIDs_ID_and_others_as_panels.txt", quote = F, col.names = F, row.names = F)