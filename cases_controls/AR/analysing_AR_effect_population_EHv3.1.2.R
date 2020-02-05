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
setwd("~/Documents/STRs/ANALYSIS/cases_controls/EHv3.1.2/")

# Load data
dedup_data = read.csv("table_STR_repeat_size_each_row_allele_EHv3.1.2_AR_CAG_simplified_dedup_050220.tsv",
                      sep = '\t',
                      header = T,
                      stringsAsFactors = F)
dim(dedup_data)
# 112108 19

# dataset 1
dedup_data_all = dedup_data
dim(dedup_data_all)
# 112108  19

# dataset 2
dedup_data_not_neuro_not_mito = dedup_data %>% 
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list), !grepl("[Mm][Ii][Tt][Oo]", panel_list)) 
dim(dedup_data_not_neuro_not_mito)
# 79345  19

# dataset 3
dedup_data_not_neuro_not_mito_not_cancer = dedup_data %>% 
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list), !grepl("[Mm][Ii][Tt][Oo]", panel_list), (programme %in% "Rare Diseases")) 
dim(dedup_data_not_neuro_not_mito_not_cancer)
# 63670  19

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
# 55376  19

# dataset1
# Check how many participants present an allele >=37
dedup_data_all %>% filter(repeat_size >= 37) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 132

dedup_data_all %>% filter(repeat_size >= 37) %>% select(participant_id) %>% pull() %>% length()
# 138

# Check how many participants present an allele >=38
dedup_data_all %>% filter(repeat_size >= 38) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 95

dedup_data_all %>% filter(repeat_size >= 38) %>% select(participant_id) %>% pull() %>% length()
# 98


# dataset2
# Check how many participants present an allele >=37
dedup_data_not_neuro_not_mito %>% filter(repeat_size >= 37) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 78

dedup_data_not_neuro_not_mito %>% filter(repeat_size >= 37) %>% select(participant_id) %>% pull() %>% length()
# 80

# Check how many participants present an allele >=38
dedup_data_not_neuro_not_mito %>% filter(repeat_size >= 38) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 52

dedup_data_not_neuro_not_mito %>% filter(repeat_size >= 38) %>% select(participant_id) %>% pull() %>% length()
# 54


# dataset3
# Check how many participants present an allele >=37
dedup_data_not_neuro_not_mito_not_cancer %>% filter(repeat_size >= 37) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 69

dedup_data_not_neuro_not_mito_not_cancer %>% filter(repeat_size >= 37) %>% select(participant_id) %>% pull() %>% length()
# 71

# Check how many participants present an allele >=38
dedup_data_not_neuro_not_mito_not_cancer %>% filter(repeat_size >= 38) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 45

dedup_data_not_neuro_not_mito_not_cancer %>% filter(repeat_size >= 38) %>% select(participant_id) %>% pull() %>% length()
# 47


# dataset 4
# Check how many participants present an allele >=37
dedup_data_pietro %>% filter(repeat_size >= 37) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 58

dedup_data_pietro %>% filter(repeat_size >= 37) %>% select(participant_id) %>% pull() %>% length()
# 62

# Check how many participants present an allele >=38
dedup_data_pietro %>% filter(repeat_size >= 38) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 44

dedup_data_pietro %>% filter(repeat_size >= 38) %>% select(participant_id) %>% pull() %>% length()
# 45
