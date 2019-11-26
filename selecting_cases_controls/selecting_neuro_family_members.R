# Objective: Select those family members (all within the same family) that have been recruited under `neuro` 
# Pilot and Main tables
# `disease_group` = neuro + `disease_subgroup` = Mitochondrial
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"
library(scatterplot3d); packageDescription("scatterplot3d", fields = "Version") # 0.3-41

# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/cases_controls/")

# Load PILOT clinical data
pilot_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                      sep = "\t",
                      header = T,
                      stringsAsFactors = F)
dim(pilot_data)
# 4974  10

# Load V7 RE MAIN clinical data
main_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_250919.tsv",
                     sep = "\t",
                     header =  T,
                     stringsAsFactors = F)
dim(main_data)
# 1056568  26

# Load the table which translates different disease subgroups
translator_table = read.csv("~/Documents/STRs/clinical_data/clinical_research_cohort/phenotyping_v140_2019-09-13_15-26-02.tsv",
                            header = T,
                            sep = "\t",
                            stringsAsFactors = F)
dim(translator_table)
# 2632  3

# 1 - Select from MAIN data all family members that have been assigned to have Neuro or Mito
# Main dataset is special, ONLY proband has enriched or populated the `specificDisease` info. Thus, we first need to have the list of familyIDs to take later all members
l_main_familyIds_neuro = main_data %>% 
  filter(disease_group %in% "Neurology and neurodevelopmental disorders" | disease_sub_group %in% "Mitochondrial disorders" | disease_sub_group %in% "mitochondrial") %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull()
length(l_main_familyIds_neuro)
# 14075


main_data_subset = main_data %>% 
  filter(rare_diseases_family_id %in% l_main_familyIds_neuro) %>%
  select(participant_id, plate_key.x, rare_diseases_family_id, participant_type, affection_status, participant_phenotypic_sex, specific_disease, disease_group, disease_sub_group, genetic_vs_reported_results, genome_build)
dim(main_data_subset)
# 655006  11

main_data_subset = unique(main_data_subset)
dim(main_data_subset)
# 35884  11

# Enrich the disease info for all members within the family - for consistency
for (i in l_main_familyIds_neuro){
  id_proband = which(main_data_subset$rare_diseases_family_id %in% i & main_data_subset$participant_type %in% "Proband")
  id_relatives = which(main_data_subset$rare_diseases_family_id %in% i & !main_data_subset$participant_type %in% "Proband")
  
  proband_specific_disease = main_data_subset$specific_disease[id_proband]
  proband_disease_group = main_data_subset$disease_group[id_proband]
  proband_disease_subgroup = main_data_subset$disease_sub_group[id_proband]
  
  if (length(id_proband) == 0){
    next
  }
  if (length(id_relatives) > 0){
    for(j in length(id_relatives)){
      main_data_subset$specific_disease[id_relatives[j]] = proband_specific_disease
      main_data_subset$disease_group[id_relatives[j]] = proband_disease_group
      main_data_subset$disease_sub_group[id_relatives[j]] = proband_disease_subgroup
    }
  }
}

dim(main_data_subset)
#

# 2 - Select from PILOT data all family members that have been assigned to have Neuro or Mito
# First we need to translate this from translator_table, to take the list of `specific_disease` we need to filter out from the pilot dataset
l_specific_disease = translator_table %>%
  filter(disease_group %in% "Neurology and neurodevelopmental disorders" | disease_subgroup %in% "Mitochondrial disorders") %>%
  select(specific_disease) %>%
  unique() %>%
  pull() 

length(l_specific_disease)  
# 25

# Pilot dataset is special, ONLY proband has enriched or populated the `specificDisease` info. Thus, we first need to have the list of familyIDs to take later all members
l_familyIds_neuro = pilot_data %>%
  filter(specificDisease %in% l_specific_disease) %>%
  select(gelFamilyId.x) %>%
  unique() %>%
  pull()
length(l_familyIds_neuro)  
# 707

pilot_data_subset = pilot_data %>%
  filter(gelFamilyId.x %in% l_familyIds_neuro) %>%
  select(gelID, plateKey, gelFamilyId.x, biological_relation_to_proband, disease_status, sex, specificDisease)
dim(pilot_data_subset)
# 1705  7

pilot_data_subset$genome = rep("GRCh37", length(pilot_data_subset$gelID))
