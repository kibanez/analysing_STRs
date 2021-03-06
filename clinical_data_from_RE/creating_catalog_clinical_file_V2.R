# Objective: create rd genomes and cancer germline info from the latest RE clinical data batch
# This time V2 - 2018/01/30
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string

library(dplyr)

setwd("~/Documents/STRs/clinical_data/clinical_data/raw/RE_clinical_data_V2/")

rd_family = read.csv("./rare_diseases_family_2020-07-07_11-24-57.tsv",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = F)
dim(rd_family)
# 20503  19

rd_participant = read.csv("./rare_diseases_participant_phen_2020-07-07_11-25-37.tsv",
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = F)
dim(rd_participant)
# 657611  23

rd_pedigree = read.csv("./rare_diseases_pedigree_member_2020-07-07_11-26-22.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)
dim(rd_pedigree)
# 58286  55

participant_info = read.csv("./participant_2020-07-07_11-23-34.tsv",
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = F)
dim(participant_info)
# 53190  97

disease_info = read.csv("./rare_diseases_participant_dise_2020-07-07_11-25-16.tsv",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = F)
dim(disease_info)
# 24497  10

path = read.csv("./sequencing_report_2020-07-07_11-23-59.tsv",
                sep = "\t",
                header = TRUE,
                stringsAsFactors = F)
dim(path)
# 32436  11

# Create genome_build
path = path %>%
  group_by(platekey) %>%
  mutate(genome_build = case_when(
    delivery_version == "V2" ~ "GRCh37",
    delivery_version == "V4" ~ "GRCh38",
    delivery_version == "unknown" ~ "unknown")) %>%
 ungroup() %>%
  as.data.frame()

path_subset = path %>% select(participant_id, platekey, type, path, genome_build)

# Let's focus on the REAL genomes we do have
all_data = path_subset %>% select(participant_id, platekey, path, type, genome_build)
dim(all_data)
# 32436  5

all_data = left_join(all_data,
                     disease_info %>% select(participant_id, disease_group, disease_sub_group, specific_disease, normalised_specific_disease),
                     by = "participant_id")
dim(all_data)
# 32727  9

all_data = left_join(all_data, 
                     participant_info %>% select(participant_id, rare_diseases_family_id, biological_relationship_to_proband, participant_type,
                                                 participant_medical_review_qc_state_code, year_of_birth, participant_phenotypic_sex, participant_karyotypic_sex, 
                                                 participant_stated_gender, programme_consent_status, programme, participant_medical_review_qc_state_code, 
                                                 participant_ethnic_category),
                     by = "participant_id")
dim(all_data)                 
# 32727  20

all_data = left_join(all_data,
                     rd_family %>% select(rare_diseases_family_id, family_group_type, family_medical_review_qc_state_code),
                     by = "rare_diseases_family_id")
dim(all_data)
# 32727  22

all_data = left_join(all_data,
                     rd_participant %>% filter(hpo_present %in% "Yes") %>% select(participant_id, hpo_term, hpo_id),
                     by = "participant_id")
dim(all_data)
# 100933  24

# population data
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      sep = ",",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59464  36

all_data = left_join(all_data,
                      popu_table %>% select(ID, best_guess_predicted_ancstry, self_reported),
                      by = c("platekey"="ID"))
dim(all_data)
# 100933  26

all_data = unique(all_data)
dim(all_data)
# 100933  26

write.table(all_data, "../../rd_genomes_all_data_130720_V2.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
