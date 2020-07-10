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
#74310    16 

rd_family = read.csv("./rare_diseases_family_2019-12-04_15-06-47.tsv",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = F)
dim(rd_family)
# 52274  11

rd_participant = read.csv("./rare_diseases_participant_phen_2019-12-04_15-07-16.tsv",
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = F)
dim(rd_participant)
# 1276264      12  

participant_info = read.csv("./participant_2019-12-04_15-09-42.tsv",
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = F)
dim(participant_info)
# 89335  59 

disease_info = read.csv("./rare_diseases_participant_dise_2019-12-04_15-08-50.tsv",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = F)
dim(disease_info)
# 40088 11 

panels_info = read.csv("./panels_applied_2019-12-04_15-10-15.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = F)
dim(panels_info)
# 202702  9 

cancer_analysis = read.csv("./cancer_analysis_2019-12-04_15-05-20.tsv",
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = F)
dim(cancer_analysis)
# 15838 77 

year_birth = read.csv("./rare_diseases_pedigree_member_2019-12-04_15-11-25.tsv",
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = F)
dim(year_birth)
# 216233  35 

platekeys = read.csv("./plated_sample_2019-12-04_15-11-53.tsv",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = F)
dim(platekeys)
# 113003  14 

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
# 111804  12 

all_data = left_join(all_data,
                     disease_info %>% select(participant_id, disease_group, disease_sub_group, specific_disease),
                     by = "participant_id")
dim(all_data)
# 113403  15 

all_data = left_join(all_data, 
                     participant_info %>% select(participant_id, participant_medical_review_qc_state_code, year_of_birth, participant_phenotypic_sex, participant_karyotypic_sex, participant_stated_gender, programme_consent_status, programme, participant_ethnic_category),
                     by = "participant_id")
dim(all_data)                 
# 113403  23

all_data = left_join(all_data,
                     rd_family %>% select(rare_diseases_family_id, family_group_type, family_medical_review_qc_state_code),
                     by = "rare_diseases_family_id")
dim(all_data)
# 113403  25

all_data = left_join(all_data,
                     panels_info %>% select(participant_id, panel_name, panel_version),
                     by = "participant_id")
dim(all_data)
# 263319  27 

all_data = left_join(all_data,
                     rd_participant %>% filter(hpo_present %in% "Yes") %>% select(participant_id, hpo_term, hpo_id),
                     by = "participant_id")
dim(all_data)
# 112318  29

all_data = left_join(all_data,
                     year_birth %>% select(participant_id, affection_status),
                     by = "participant_id")
dim(all_data)
# 1124633  30

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
# 1124633  32

write.table(all_data, "../../rd_genomes_all_data_100720_V8.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
