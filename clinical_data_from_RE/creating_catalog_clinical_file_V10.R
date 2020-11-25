# Objective: create rd genomes and cancer germline info from the latest RE clinical data batch
# Update on 07/10/20: including some fields from `clinic_sample` and `gmc_exit_questionnaire` tables
# Update on 24/11/20: include secondary information from HPC
# This time V10 - 2020/09/03
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string

library(dplyr)

setwd("~/Documents/STRs/clinical_data/clinical_data/raw/RE_clinical_data_V10")

rd_analysis = read.csv("./rare_disease_analysis_2020-09-08_09-28-42.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)

dim(rd_analysis)
# 74271  16

rd_family = read.csv("./rare_diseases_family_2020-09-08_09-37-49.tsv",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = F)
dim(rd_family)
# 52981  11

rd_participant = read.csv("./rare_diseases_participant_phen_2020-09-08_09-33-20.tsv",
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = F)
dim(rd_participant)
# 1274933  14

participant_info = read.csv("./participant_2020-09-08_09-31-01.tsv",
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = F)
dim(participant_info)
# 89139  59

disease_info = read.csv("./rare_diseases_participant_dise_2020-09-08_09-35-50.tsv",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = F)
dim(disease_info)
# 39789  11

panels_info = read.csv("./panels_applied_2020-09-08_09-31-27.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)
dim(panels_info)
# 207966  9

year_birth = read.csv("./rare_diseases_pedigree_member_2020-09-08_09-36-18.tsv",
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = F)
dim(year_birth)
# 217473  35

path = read.csv("./genome_file_paths_and_types_2020-09-08_09-36-43.tsv",
                sep = "\t",
                header = TRUE,
                stringsAsFactors = F)
dim(path)
# 518456  12

clinic_sample = read.csv("./clinic_sample_2020-10-07_12-09-30.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(clinic_sample)
# 190647  55

gmc_exit = read.csv("./gmc_exit_questionnaire_2020-10-07_13-07-38.tsv",
                    sep = "\t",
                    stringsAsFactors = F, 
                    header = T)
dim(gmc_exit)
# 28757  23

path_subset = path %>% filter(file_sub_type %in% "BAM") %>% select(participant_id, platekey, type, file_path)
rm(path)

# Let's focus on the REAL genomes we do have
all_data = path_subset %>% select(participant_id, platekey, file_path, type)
dim(all_data)
# 107540  4
rm(path_subset)

all_data = left_join(all_data,
                     rd_analysis %>% select(participant_id, rare_diseases_family_id, plate_key, biological_relationship_to_proband, participant_type, normalised_specific_disease, genome_build, genetic_vs_reported_results, participant_ethnic_category),
                     by = "participant_id")
dim(all_data)
# 112021  12
rm(rd_analysis)

all_data = left_join(all_data,
                     disease_info %>% select(participant_id, disease_group, disease_sub_group, specific_disease),
                     by = "participant_id")
dim(all_data)
# 113538  15
rm(disease_info)

all_data = left_join(all_data, 
                     participant_info %>% select(participant_id, registered_at_gmc_trust, participant_medical_review_qc_state_code, year_of_birth, participant_phenotypic_sex, participant_karyotypic_sex, participant_stated_gender, programme_consent_status, programme),
                     by = "participant_id")
dim(all_data)                 
# 113538  23
rm(participant_info)

all_data = left_join(all_data,
                     rd_family %>% select(rare_diseases_family_id, family_group_type, family_medical_review_qc_state_code),
                     by = "rare_diseases_family_id")
dim(all_data)
# 113538  25
rm(rd_family)

all_data = left_join(all_data,
                     panels_info %>% select(participant_id, panel_name, panel_version),
                     by = "participant_id")
dim(all_data)
# 272617  27
rm(panels_info)

all_data = left_join(all_data,
                     rd_participant %>% filter(hpo_present %in% "Yes") %>% select(participant_id, hpo_term, hpo_id),
                     by = "participant_id")
dim(all_data)
# 1181424  29
rm(rd_participant)

all_data = left_join(all_data,
                     year_birth %>% select(participant_id, affection_status),
                     by = "participant_id")
dim(all_data)
# 1184730  30
rm(year_birth)

all_data = left_join(all_data,
                     clinic_sample %>% select(participant_id, clinic_sample_collected_at_gmc, clinic_sample_collected_at_gmc_trust),
                     by = "participant_id")
dim(all_data)
# 2399541  32
rm(clinic_sample)    

all_data = left_join(all_data,
                     gmc_exit %>% select(participant_id, case_solved_family),
                     by = "participant_id")
dim(all_data)
# 3474081  33
rm(gmc_exit)

# New RESCTY and POSTDIST info for each PID
hpc = read.csv("./hes_op_2020-11-24_09-11-57_pid_postdist_rescty.tsv",
               sep = "\t",
               stringsAsFactors = F,
               header = T)
dim(hpc)
# 5628380  3

hpc = unique(hpc)
dim(hpc)
# 145770  3

all_data = left_join(all_data,
                     hpc,
                     by = "participant_id")
dim(all_data)
# 6544650  35
rm(hpc)

all_data = unique(all_data)
dim(all_data)
#

# population data - let's enrich with merged (batch1 and batch2) population info
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/popu_merged_batch1_batch2_79849_genomes.tsv",
                      sep = "\t",
                      stringsAsFactors = F,
                      header = F)
dim(popu_table)
# 79849  2
colnames(popu_table) = c("platekey", "superpopu")

all_data = left_join(all_data,
                      popu_table,
                      by = "platekey")
dim(all_data)
# 3474081 33

write.table(all_data, "../../rd_genomes_all_data_251120_V10.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
