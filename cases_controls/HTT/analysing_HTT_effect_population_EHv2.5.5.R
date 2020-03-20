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
setwd("~/Documents/STRs/ANALYSIS/cases_controls/EHv2.5.5/")

# Load data
dedup_data = read.csv("table_STR_repeat_size_each_row_allele_EHv2.5.5_HTT_CAG_simplified_dedup_040220.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(dedup_data)
# 132102  19


# New thing I want to retrieve from here in terms of population for this 60K unique participants
# And also plot the distribution of gender and age within this subset
dedup_data = dedup_data %>%
  group_by(participant_id) %>%
  mutate(age = 2020 - year_of_birth) %>%
  ungroup() %>%
  as.data.frame()
  
age_distribution = ggplot(dedup_data,
                          aes(x = age, fill = participant_phenotypic_sex)) + 
  geom_density(alpha=0.25) +
  xlab("age of participant") +
  ylab("density") +
  scale_fill_discrete(name = 'gender')

png("HTT_distribution_gender_age_66051_participants.png")
print(age_distribution)
dev.off()

popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/matthias_work_main/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      sep = ",",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59464  36

popu_subset = popu_table %>%
  filter(ID %in% unique(dedup_data$platekey))
dim(popu_subset)
# 50626  36



# Definition of different sub-datasets in `dedup_data`

# dataset 1
dedup_data_all = dedup_data
dim(dedup_data_all)
# 132102  19

# dataset 2
dedup_data_not_neuro_not_mito = dedup_data %>% 
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list), !grepl("[Mm][Ii][Tt][Oo]", panel_list)) 
dim(dedup_data_not_neuro_not_mito)
# 86894  19

# dataset 3
dedup_data_not_neuro_not_mito_not_cancer = dedup_data %>% 
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list), !grepl("[Mm][Ii][Tt][Oo]", panel_list), (programme %in% "Rare Diseases")) 
dim(dedup_data_not_neuro_not_mito_not_cancer)
# 84572  19

dedup_only_cancer = dedup_data %>%
  filter(programme %in% "Cancer")
dim(dedup_only_cancer)
# 2322  19

# dataset 4 - only probands
dedup_only_probands = dedup_data %>% 
  filter((biological_relationship_to_proband %in% "N/A" & programme %in% "Rare Diseases") | programme %in% "Cancer")
dim(dedup_only_probands)
# 62894  19

# dataset 5 - probands minus cancer
dedup_only_probands_minus_cancer = dedup_data %>% 
  filter(biological_relationship_to_proband %in% "N/A" & programme %in% "Rare Diseases")
dim(dedup_only_probands_minus_cancer)
# 60572  19

# dataset 6 - probands minus cancer minus neurology
dedup_only_probands_minus_cancer_minus_neuro = dedup_data %>% 
  filter(biological_relationship_to_proband %in% "N/A" & programme %in% "Rare Diseases" & !grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list) & !grepl("[Mm][Ii][Tt][Oo]", panel_list))
dim(dedup_only_probands_minus_cancer_minus_neuro)
# 33084  19

# dataset 7 - probands minus neuro
dedup_only_probands_minus_neuro = dedup_data %>% 
  filter((biological_relationship_to_proband %in% "N/A" & programme %in% "Rare Diseases" & !grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group_list) & !grepl("[Mm][Ii][Tt][Oo]", panel_list)) | programme %in% "Cancer")
dim(dedup_only_probands_minus_neuro)
# 35406  19


# Let's compute numbers not taking into account RELATEDNESS
# All participants
# patho (>=40)
dedup_data_all %>% filter(repeat_size >= 40) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 30
dedup_data_all %>% filter(repeat_size >= 40) %>% select(participant_id) %>% pull() %>% length()
# 30

# reduced penetrance (36-39)
dedup_data_all %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 134
dedup_data_all %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% pull() %>% length()
# 135

# intermediate (27-35)
dedup_data_all %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 4197
dedup_data_all %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% pull() %>% length()
# 4306


# ONLY probands
# patho (>=40)
dedup_only_probands %>% filter(repeat_size >= 40) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 19
dedup_only_probands %>% filter(repeat_size >= 40) %>% select(participant_id) %>% pull() %>% length()
# 19

# reduced penetrance (36-39)
dedup_only_probands %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 58
dedup_only_probands %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% pull() %>% length()
# 59

# intermediate (27-35)
dedup_only_probands %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 1963
dedup_only_probands %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% pull() %>% length()
# 2030


# Probands minus cancer
# patho (>=40)
dedup_only_probands_minus_cancer %>% filter(repeat_size >= 40) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 18
dedup_only_probands_minus_cancer %>% filter(repeat_size >= 40) %>% select(participant_id) %>% pull() %>% length()
# 18

# reduced penetrance (36-39)
dedup_only_probands_minus_cancer %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 54
dedup_only_probands_minus_cancer %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% pull() %>% length()
# 55

# intermediate (27-35)
dedup_only_probands_minus_cancer %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 1897
dedup_only_probands_minus_cancer %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% pull() %>% length()
# 1962


# Probands minus neuro
# patho (>=40)
dedup_only_probands_minus_neuro %>% filter(repeat_size >= 40) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 7
dedup_only_probands_minus_neuro %>% filter(repeat_size >= 40) %>% select(participant_id) %>% pull() %>% length()
# 7

# reduced penetrance (36-39)
dedup_only_probands_minus_neuro %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 35
dedup_only_probands_minus_neuro %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% pull() %>% length()
# 35

# intermediate (27-35)
dedup_only_probands_minus_neuro %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 1156
dedup_only_probands_minus_neuro %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% pull() %>% length()
# 1186

# Probands minus neuro and cancer
# patho (>=40)
dedup_only_probands_minus_cancer_minus_neuro %>% filter(repeat_size >= 40) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 6
dedup_only_probands_minus_cancer_minus_neuro %>% filter(repeat_size >= 40) %>% select(participant_id) %>% pull() %>% length()
# 6

# reduced penetrance (36-39)
dedup_only_probands_minus_cancer_minus_neuro %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 31
dedup_only_probands_minus_cancer_minus_neuro %>% filter(repeat_size >= 36 & repeat_size <=39) %>% select(participant_id) %>% pull() %>% length()
# 31

# intermediate (27-35)
dedup_only_probands_minus_cancer_minus_neuro %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 1090
dedup_only_probands_minus_cancer_minus_neuro %>% filter(repeat_size >= 27 & repeat_size <=35) %>% select(participant_id) %>% pull() %>% length()
# 1118


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
# 285078  3

# Let's visualise them with violin plots

ggplot(merge_all, aes(x = dataset_name, y=repeat_size, fill = dataset_name)) +
  geom_violin() +
  xlab("Defined datasets") + 
  ylab("Repeat sizes (repeat units)") 


# Taking repeat-size >= 40
violin_beyond40 = ggplot(merge_all %>% filter(repeat_size >= 40),
       aes(x = dataset_name, y=repeat_size, fill = dataset_name)) +
  geom_violin() +
  xlab("Defined datasets, repeat-sizes >= 40") + 
  ylab("Repeat sizes (repeat units)") 

# Taking repeat-size >= 40 AND age beyond 35
violin_beyond40_age35 = ggplot(merge_all %>% filter(repeat_size >= 40, year_of_birth <= 1985),
                         aes(x = dataset_name, y=repeat_size, fill = dataset_name)) +
  geom_violin() +
  xlab("Defined datasets, repeat-sizes >= 40 and year <= 1985") + 
  ylab("Repeat sizes (repeat units)") 

# 36<= repeat-size <=39
violin_bet36and39 = ggplot(merge_all %>% filter(repeat_size >= 36 | repeat_size <= 39),
       aes(x = dataset_name, y=repeat_size, fill = dataset_name)) +
  geom_violin() +
  xlab("Defined datasets, 36<=repeat-sizes<=39") + 
  ylab("Repeat sizes (repeat units)") 

# 27<= repeat-size <=35
violin_bet27and35 = ggplot(merge_all %>% filter(repeat_size >= 27 | repeat_size <= 35),
       aes(x = dataset_name, y=repeat_size, fill = dataset_name)) +
  geom_violin() +
  xlab("Defined datasets, 27<=repeat-sizes<=35") + 
  ylab("Repeat sizes (repeat units)") 


# save the violin plots
setwd("./analyses/HTT")

png("distribution_across_datasets_beyond40.png")
print(violin_beyond40)
dev.off()

png("distribution_across_datasets_beyond40_beyond35.png")
print(violin_beyond40_age35)
dev.off()

png("distribution_across_datasets_bet27and35.png")
print(violin_bet27and35)
dev.off()

png("distribution_across_datasets_bet36and39.png")
print(violin_bet36and39)
dev.off()

# Number of genomes/participants having repeat-size>=40 , so on and so forth


