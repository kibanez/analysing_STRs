# Objective: create rd genomes and cancer germline info from the latest RE clinical data batch
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string


library(dplyr)


setwd("~/Downloads/")

rd_analysis = read.csv("./rare_disease_analysis_2019-09-17_14-53-35.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)

dim(rd_analysis)
#74981 16

rd_family = read.csv("./rare_diseases_family_2019-09-25_10-21-58.tsv",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = F)
dim(rd_family)
#49396  11

rd_participant = read.csv("./rare_diseases_participant_phen_2019-09-17_14-56-19.tsv",
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = F)
dim(rd_participant)
#1288624 12


participant_info = read.csv("./participant_2019-09-17_14-54-59.tsv",
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = F)
dim(participant_info)
#90643 59

disease_info = read.csv("./rare_diseases_participant_dise_2019-09-18_11-06-19.tsv",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = F)
dim(disease_info)
#40611 10

panels_info = read.csv("./panels_applied_2019-09-17_14-55-37.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)
dim(panels_info)
#183861 8


cancer_analysis = read.csv("./cancer_analysis_2019-09-17_14-54-28.tsv",
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = F)
dim(cancer_analysis)
#12594  77

year_birth = read.csv("./rare_diseases_pedigree_member_2019-09-17_16-45-42.tsv",
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = F)
dim(year_birth)
#206585   35

platekeys = read.csv("./plated_sample_2019-09-20_11-35-13.tsv",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = F)
dim(platekeys)
#113696 13

all_data = platekeys %>% select(participant_id, plate_key)
dim(all_data)
#113696 2

all_data = left_join(all_data,
                     rd_analysis %>% select(participant_id, rare_diseases_family_id, plate_key, biological_relationship_to_proband, participant_type, normalised_specific_disease, genome_build, genetic_vs_reported_results),
                     by = "participant_id")
dim(all_data)
#115988  9

all_data = left_join(all_data,
                     disease_info %>% select(participant_id, disease_group, disease_sub_group, specific_disease),
                     by = "participant_id")
dim(all_data)
#117542 12


all_data = left_join(all_data, 
                     participant_info %>% select(participant_id, participant_medical_review_qc_state_code, year_of_birth, participant_phenotypic_sex, participant_karyotypic_sex, participant_stated_gender, programme_consent_status, programme),
                     by = "participant_id")
dim(all_data)                 
#117542 19


all_data = left_join(all_data,
                     rd_family %>% select(rare_diseases_family_id, family_group_type, family_medical_review_qc_state_code),
                     by = "rare_diseases_family_id")
dim(all_data)
#117542 21


all_data = left_join(all_data,
                     panels_info %>% select(participant_id, panel_name, panel_version),
                     by = "participant_id")
dim(all_data)
#253170 23

all_data = left_join(all_data,
                     rd_participant %>% filter(hpo_present %in% "Yes") %>% select(participant_id, hpo_term, hpo_id),
                     by = "participant_id")
dim(all_data)
#1055635 25

all_data = left_join(all_data,
                     year_birth %>% select(participant_id, affection_status),
                     by = "participant_id")
dim(all_data)
#1056568 26

write.table(all_data, "rd_genomes_all_data_250919.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

