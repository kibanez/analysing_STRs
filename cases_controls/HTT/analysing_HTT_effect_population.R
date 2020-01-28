# Objective: is to see the distribution of repeat-expansions on HTT across 100K
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
dedup_data = read.csv("table_STR_repeat_size_each_row_allele_EHv3.1.2_HTT_CAG_simplified_dedup_270120.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(dedup_data)
# 130126  19

# Definition of different sub-datasets in `dedup_data`

# dataset 1
dedup_data_all = dedup_data
dim(dedup_data_all)
# 130126  19

# dataset 2
dedup_data_not_neuro_not_mito = dedup_data %>% 
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list), !grepl("[Mm][Ii][Tt][Oo]", panel_list)) 
dim(dedup_data_not_neuro_not_mito)
# 86540  19

# dataset 3
dedup_data_not_neuro_not_mito_not_cancer = dedup_data %>% 
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list), !grepl("[Mm][Ii][Tt][Oo]", panel_list), (programme %in% "Rare Diseases")) 
dim(dedup_data_not_neuro_not_mito_not_cancer)
# 81602  19

# For each dataset, we ONLY want to consider `UNRELATED`genomes
# STRATEGY 1 - take as UNRELATED genomes coming from `probands`

# dataset 1 (cancer programme I assume they are unrelated)
dedup_data_all_unrelated1 = dedup_data_all %>%
  filter((biological_relationship_to_proband %in% "N/A" & programme %in% "Rare Diseases") | programme %in% "Cancer")
dim(dedup_data_all_unrelated1)
# 63526  19

length(unique(dedup_data_all_unrelated1$participant_id))
# 31763
length(unique(dedup_data_all_unrelated1$platekey))
# 31763
length(unique(dedup_data_all_unrelated1$rare_diseases_family_id))
# 29295
length(which(is.na(dedup_data_all_unrelated1$rare_diseases_family_id)))
# 4938

# dataset 2 
dedup_data_not_neuro_not_mito_unrelated1 = dedup_data_not_neuro_not_mito %>%
  filter((biological_relationship_to_proband %in% "N/A" & programme %in% "Rare Diseases") | programme %in% "Cancer")
dim(dedup_data_not_neuro_not_mito_unrelated1)
# 36994  19

length(unique(dedup_data_not_neuro_not_mito_unrelated1$participant_id))
# 18497
length(unique(dedup_data_not_neuro_not_mito_unrelated1$platekey))
# 18497
length(unique(dedup_data_not_neuro_not_mito_unrelated1$rare_diseases_family_id))
# 16029
length(which(is.na(dedup_data_not_neuro_not_mito_unrelated1$rare_diseases_family_id)))
# 4938

# dataset 3
dedup_data_not_neuro_not_mito_not_cancer_unrelated1 = dedup_data_not_neuro_not_mito_not_cancer %>%
  filter(biological_relationship_to_proband %in% "N/A" & programme %in% "Rare Diseases")
dim(dedup_data_not_neuro_not_mito_not_cancer_unrelated1)
# 32056  19

length(unique(dedup_data_not_neuro_not_mito_not_cancer_unrelated1$participant_id))
# 16028
length(unique(dedup_data_not_neuro_not_mito_not_cancer_unrelated1$platekey))
# 16028
length(unique(dedup_data_not_neuro_not_mito_not_cancer_unrelated1$rare_diseases_family_id))
# 16028


# STRATEGY 2 - take as UNRELEATED genomes that are `Father` OR `Mother` (take both if they are available, we could have ~1% of relatedness....)

# dataset 1
dedup_data_all_unrelated2 = dedup_data_all %>%
  filter((biological_relationship_to_proband %in% c("Father", "Mother") & programme %in% "Rare Diseases") | programme %in% "Cancer")
dim(dedup_data_all_unrelated2)  
# 62128  19

length(unique(dedup_data_all_unrelated2$participant_id))
# 31064
length(unique(dedup_data_all_unrelated2$platekey))
# 31064
length(unique(dedup_data_all_unrelated2$rare_diseases_family_id))
# 17289
length(which(is.na(dedup_data_all_unrelated2$rare_diseases_family_id)))
# 4938

# dataset 2
dedup_data_not_neuro_not_mito_unrelated2 = dedup_data_not_neuro_not_mito %>%
  filter((biological_relationship_to_proband %in% c("Father", "Mother") & programme %in% "Rare Diseases") | programme %in% "Cancer")
dim(dedup_data_not_neuro_not_mito_unrelated2)  
# 47656  19

length(unique(dedup_data_not_neuro_not_mito_unrelated2$participant_id))
# 23828
length(unique(dedup_data_not_neuro_not_mito_unrelated2$platekey))
# 23828
length(unique(dedup_data_not_neuro_not_mito_unrelated2$rare_diseases_family_id))
# 13065
length(which(is.na(dedup_data_not_neuro_not_mito_unrelated2$rare_diseases_family_id)))
# 4938

# dataset 3
dedup_data_not_neuro_not_mito_not_cancer_unrelated2 = dedup_data_not_neuro_not_mito_not_cancer %>%
  filter(biological_relationship_to_proband %in% c("Father", "Mother") & programme %in% "Rare Diseases")
dim(dedup_data_not_neuro_not_mito_not_cancer_unrelated2)  
# 42718  19

length(unique(dedup_data_not_neuro_not_mito_not_cancer_unrelated2$participant_id))
# 21359
length(unique(dedup_data_not_neuro_not_mito_not_cancer_unrelated2$platekey))
# 21359
length(unique(dedup_data_not_neuro_not_mito_not_cancer_unrelated2$rare_diseases_family_id))
# 13064


# Once we have defined all 3 datasets to go through our analysis, let's see the differences in the distribution of str-repeats across them
# Let's merge all them, just taking the repeat-size and the group name
dataset1_stra1 = dedup_data_all_unrelated1 %>%
  mutate(dataset_name = "dataset1_stra1") %>%
  select(repeat_size, year_of_birth, dataset_name)

dataset2_stra1 = dedup_data_not_neuro_not_mito_unrelated1 %>%
  mutate(dataset_name = "dataset2_stra1") %>%
  select(repeat_size, year_of_birth, dataset_name)

dataset3_stra1 = dedup_data_not_neuro_not_mito_not_cancer_unrelated1 %>%
  mutate(dataset_name = "dataset3_stra1") %>%
  select(repeat_size, year_of_birth, dataset_name)

dataset1_stra2 = dedup_data_all_unrelated2 %>%
  mutate(dataset_name = "dataset1_stra2") %>%
  select(repeat_size, year_of_birth, dataset_name)

dataset2_stra2 = dedup_data_not_neuro_not_mito_unrelated2 %>%
  mutate(dataset_name = "dataset2_stra2") %>%
  select(repeat_size, year_of_birth, dataset_name)

dataset3_stra2 = dedup_data_not_neuro_not_mito_not_cancer_unrelated2 %>%
  mutate(dataset_name = "dataset3_stra2") %>%
  select(repeat_size, year_of_birth, dataset_name)

merge_all = rbind(dataset1_stra1,
                  dataset1_stra2,
                  dataset2_stra1,
                  dataset2_stra2,
                  dataset3_stra1,
                  dataset3_stra2)
dim(merge_all)
# 285078  2

# Let's visualise them with violin plots

ggplot(merge_all, aes(x = dataset_name, y=repeat_size, fill = dataset_name)) +
  geom_violin() +
  xlab("Defined datasets") + 
  ylab("Repeat sizes (repeat units)") 


# Taking repeat-size >= 40
