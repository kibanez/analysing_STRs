# Objective: is to see the distribution of repeat-expansions on AR across 100K
# 1) Definining 3 different threshold cutoffs
# 2) Definining 3 different sub-cohorts or sub-datasets within 85K
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/EHv2.5.5/")

# Load data
dedup_data = read.csv("table_STR_repeat_size_each_row_allele_EHv2.5.5_AR_CAG_simplified_dedup_050220.tsv",
                      sep = '\t',
                      header = T,
                      stringsAsFactors = F)
dim(dedup_data)
# 115772 19

# dataset 1
dedup_data_all = dedup_data
dim(dedup_data_all)
# 115772  19

# dataset 2
dedup_data_not_neuro_not_mito = dedup_data %>% 
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list), !grepl("[Mm][Ii][Tt][Oo]", panel_list)) 
dim(dedup_data_not_neuro_not_mito)
# 81790  19

# dataset 3
dedup_data_not_neuro_not_mito_not_cancer = dedup_data %>% 
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list), !grepl("[Mm][Ii][Tt][Oo]", panel_list), (programme %in% "Rare Diseases")) 
dim(dedup_data_not_neuro_not_mito_not_cancer)
# 66003  19

dedup_only_cancer = dedup_data %>%
  filter(programme %in% "Cancer")
dim(dedup_only_cancer)
# 15787  19

# dataset 4 - female probands OR male probands NOT cancer
dedup_data_pietro = dedup_data %>%
  filter((biological_relationship_to_proband %in% "N/A" & participant_phenotypic_sex %in% "Female") |
         (is.na(biological_relationship_to_proband) & participant_phenotypic_sex %in% "Female") |
         (biological_relationship_to_proband %in% "N/A" & participant_phenotypic_sex %in% "Male" & programme %in% "Rare Diseases") |
         (is.na(biological_relationship_to_proband) & participant_phenotypic_sex %in% "Male" & programme %in% "Rare Diseases"))
dim(dedup_data_pietro)
# 56951  19


# Check how many participants present an allele >=37
dedup_data_pietro %>% filter(repeat_size >= 37) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 44

dedup_data_pietro %>% filter(repeat_size >= 37) %>% select(participant_id) %>% pull() %>% length()
# 46

# Check how many participants present an allele >=38
dedup_data_pietro %>% filter(repeat_size >= 38) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 32

dedup_data_pietro %>% filter(repeat_size >= 38) %>% select(participant_id) %>% pull() %>% length()
# 33
