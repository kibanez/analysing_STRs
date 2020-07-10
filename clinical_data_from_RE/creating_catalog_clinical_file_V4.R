# Objective: create rd genomes and cancer germline info from the latest RE clinical data batch
# This time V4 - 2018/07/31
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string

library(dplyr)

setwd("~/Documents/STRs/clinical_data/clinical_data/raw/RE_clinical_data_V4/")

rd_analysis = read.csv("./rare_disease_analysis_2020-07-07_10-55-22.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)

dim(rd_analysis)
# 31972  12

rd_family = read.csv("./rare_diseases_family_2020-07-07_10-56-59.tsv",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = F)
dim(rd_family)
# 28920  11

rd_participant = read.csv("./rare_diseases_participant_phen_2020-07-07_10-58-13.tsv",
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = F)
dim(rd_participant)
# 853440  11

participant_info = read.csv("./participant_2020-07-07_10-56-10.tsv",
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = F)
dim(participant_info)
# 71331  55

disease_info = read.csv("./rare_diseases_participant_dise_2020-07-07_10-57-49.tsv",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = F)
dim(disease_info)
# 31938  10

panels_info = read.csv("./panels_applied_2020-07-07_10-57-28.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)
dim(panels_info)
# 57644  4

cancer_analysis = read.csv("./cancer_analysis_2020-07-07_10-55-01.tsv",
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = F)
dim(cancer_analysis)
# 20144  31

path = read.csv("./genome_file_paths_and_types_2020-07-07_10-55-46.tsv",
                sep = "\t",
                header = TRUE,
                stringsAsFactors = F)
dim(path)
# 297181  10

path_subset = path %>% filter(file_sub_type %in% "BAM") %>% select(participant_id, platekey, type, file_path)

# Let's focus on the REAL genomes we do have
all_data = path_subset %>% select(participant_id, platekey, file_path, type)
dim(all_data)
# 55472  4

all_data = left_join(all_data,
                     rd_analysis %>% select(participant_id, rare_diseases_family_id, platekey, biological_relationship_to_proband, participant_type, normalised_specific_disease, genome_build, participant_ethnic_category),
                     by = "participant_id")
dim(all_data)
# 80547  11

all_data = left_join(all_data,
                     disease_info %>% select(participant_id, disease_group, disease_sub_group, specific_disease),
                     by = "participant_id")
dim(all_data)
# 81554  14

all_data = left_join(all_data, 
                     participant_info %>% select(participant_id, participant_medical_review_qc_state_code, year_of_birth, participant_phenotypic_sex, participant_karyotypic_sex, participant_stated_gender, programme_consent_status, programme, participant_medical_review_qc_state_code),
                     by = "participant_id")
dim(all_data)                 
# 57040  22

all_data = left_join(all_data,
                     rd_family %>% select(rare_diseases_family_id, family_group_type, family_medical_review_qc_state_code),
                     by = "rare_diseases_family_id")
dim(all_data)
# 57040  24

all_data = left_join(all_data,
                     panels_info %>% select(participant_id, panel_name, panel_version),
                     by = "participant_id")
dim(all_data)
# 95981  26

all_data = left_join(all_data,
                     rd_participant %>% filter(hpo_present %in% "Yes") %>% select(participant_id, hpo_term, hpo_id),
                     by = "participant_id")
dim(all_data)
# 370296  28

# population data
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      sep = ",",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59464  36

all_data = left_join(all_data,
                      popu_table %>% select(ID, best_guess_predicted_ancstry, self_reported),
                      by = c("platekey.x"="ID"))
dim(all_data)
# 370296  30

write.table(all_data, "../../rd_genomes_all_data_100720_V4.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
