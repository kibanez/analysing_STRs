# Objective: create rd genomes and cancer germline info from the latest RE clinical data batch
# This time V8 - 28/11/2019
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string

library(dplyr)

setwd("~/Documents/STRs/clinical_data/clinical_data/raw/RE_clinical_data_V8/")

rd_analysis = read.csv("./rare_disease_analysis_2019-12-04_15-03-51.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)

dim(rd_analysis)
#74981 16 (V7)
#74310    16 (V8, now)

rd_family = read.csv("./rare_diseases_family_2019-12-04_15-06-47.tsv",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = F)
dim(rd_family)
#49396  11 (v7)
# 52274  11 (V8)

rd_participant = read.csv("./rare_diseases_participant_phen_2019-12-04_15-07-16.tsv",
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = F)
dim(rd_participant)
#1288624 12 (V7)
#1276264      12  (V8)


participant_info = read.csv("./participant_2019-12-04_15-09-42.tsv",
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = F)
dim(participant_info)
#90643 59 (V7)
#89335  59 (V8) 

disease_info = read.csv("./rare_diseases_participant_dise_2019-12-04_15-08-50.tsv",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = F)
dim(disease_info)
#40611 10 (V7)
# 40088 11 (V8)

panels_info = read.csv("./panels_applied_2019-12-04_15-10-15.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)
dim(panels_info)
#183861 8 (V7)
# 202702  9 (V8)


cancer_analysis = read.csv("./cancer_analysis_2019-12-04_15-05-20.tsv",
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = F)
dim(cancer_analysis)
#12594  77 (V7)
# 15838 77 (V8)

year_birth = read.csv("./rare_diseases_pedigree_member_2019-12-04_15-11-25.tsv",
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = F)
dim(year_birth)
#206585   35 (V7)
# 216233  35 (V8)

platekeys = read.csv("./plated_sample_2019-12-04_15-11-53.tsv",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = F)
dim(platekeys)
#113696 13 (V7)
# 113003  14 (V8)

path = read.csv("./genome_file_paths_and_types_2019-12-04_15-13-29.tsv",
                sep = "\t",
                header = TRUE,
                stringsAsFactors = F)
dim(path)
# 500443  11

path_subset = path %>% filter(file_sub_type %in% "BAM") %>% select(participant_id, platekey, type, file_path)

# Let's focus on the REAL genomes we do have
all_data = path_subset %>% select(participant_id, platekey, file_path, type)
dim(all_data)
# 107623  4

all_data = left_join(all_data,
                     rd_analysis %>% select(participant_id, rare_diseases_family_id, plate_key, biological_relationship_to_proband, participant_type, normalised_specific_disease, genome_build, genetic_vs_reported_results, participant_ethnic_category),
                     by = "participant_id")
dim(all_data)
#115988  9 (V7)
# 111804  11 (V8)

all_data = left_join(all_data,
                     disease_info %>% select(participant_id, disease_group, disease_sub_group, specific_disease),
                     by = "participant_id")
dim(all_data)
# 117542 12 (V7)
# 113403  14 (V8)

all_data = left_join(all_data, 
                     participant_info %>% select(participant_id, participant_medical_review_qc_state_code, year_of_birth, participant_phenotypic_sex, participant_karyotypic_sex, participant_stated_gender, programme_consent_status, programme),
                     by = "participant_id")
dim(all_data)                 
# 117542 19 (V7)
# 113403  21 (V8)


all_data = left_join(all_data,
                     rd_family %>% select(rare_diseases_family_id, family_group_type, family_medical_review_qc_state_code),
                     by = "rare_diseases_family_id")
dim(all_data)
# 117542 21 (V7)
# 113403  23 (V8)


all_data = left_join(all_data,
                     panels_info %>% select(participant_id, panel_name, panel_version),
                     by = "participant_id")
dim(all_data)
# 253170 23 (v7)
# 263319  25 (V8)

all_data = left_join(all_data,
                     rd_participant %>% filter(hpo_present %in% "Yes") %>% select(participant_id, hpo_term, hpo_id),
                     by = "participant_id")
dim(all_data)
# 1055635 25 (V7)
# 112318  27 (V8)

all_data = left_join(all_data,
                     year_birth %>% select(participant_id, affection_status),
                     by = "participant_id")
dim(all_data)
# 1056568 26 (V7)
# 1124633  28 (V8)

# population data
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/matthias_work_main/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      sep = ",",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59464  36

all_data = left_join(all_data,
                      popu_table %>% select(ID, best_guess_predicted_ancstry),
                      by = c("platekey"="ID"))
dim(all_data)
# 1124633  29



write.table(all_data, "../../rd_genomes_all_data_230320.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
