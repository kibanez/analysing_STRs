# Objective: merge Pilot and Main programmes' clinical data
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string

library(dplyr)

# Set working directory
setwd("~/Documents/STRs/clinical_data/clinical_data/")

# Load latest Pilot data (frozen)
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           header = TRUE)
dim(pilot_clin_data)
# 4974  10

# Let´s put all panel names into 1 single string splitted by ','
list_panels_pilot = pilot_clin_data %>% group_by(gelID) %>% summarise(panel_list = toString(unique(specificDisease))) %>% ungroup() %>% as.data.frame()
dim(list_panels_pilot)
# 4833  2

pilot_clin_data = left_join(pilot_clin_data,
                            list_panels_pilot,
                            by = "gelID")
# Remove specificDisease
pilot_clin_data = pilot_clin_data[,-8]
pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 4833  10

# Let's enrich with popu data
pilot_popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            header = T,
                            sep = ",")
dim(pilot_popu_table)
# 4821  44

pilot_clin_data = left_join(pilot_clin_data,
                            pilot_popu_table %>% select(ID, bestGUESS_sub_pop, bestGUESS_super_pop, self_reported),
                            by = c("plateKey"="ID"))

pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 4834  13

# Enrich pilot clinical data with disease group and disease subgroup
pheno_table = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/phenotyping_v140_2019-09-13_15-26-02.tsv",
                       stringsAsFactors = F,
                       header = T,
                       sep = "\t")
pheno_table = unique(pheno_table)
dim(pheno_table)
# 106  3

pilot_clin_data = left_join(pilot_clin_data,
                            pheno_table,
                            by = c("panel_list" = "specific_disease"))
pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 4834  15

# Load latests Main clinical data release (created from our R scripts)
# On 14/04/2021 I realised there were many PIDs in RE V7 not included in V10 nor V11
clin_data_V7 = read.table("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_100720_V7.tsv",
                       sep = "\t",
                       stringsAsFactors = FALSE, 
                       header = TRUE)
dim(clin_data_V7)  
# 1021801  31 

l_unique_PIDs_V7 = unique(clin_data_V7$participant_id)
length(l_unique_PIDs_V7)
# 86682

clin_data = read.table("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_251120_V10.tsv",
                       sep = "\t",
                       stringsAsFactors = FALSE, 
                       header = TRUE)
dim(clin_data)  
# 2096500 36

clin_data2 = read.table("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_301220_V11.tsv",
                       sep = "\t",
                       stringsAsFactors = FALSE, 
                       header = TRUE)
dim(clin_data2)  
# 2119961  36

# There are 164 PID in V10 and not in V11, and 111 in V11 and not in V10
clin_data = rbind(clin_data,
                  clin_data2)
clin_data = unique(clin_data)
dim(clin_data)
# 2440099  36

# Which PIDs are in V7 and not up to V11?
l_missing_fromV7 = setdiff(l_unique_PIDs_V7, unique(clin_data$participant_id))
length(l_missing_fromV7)
# 1373

clin_data_V7 = clin_data_V7 %>% 
  filter(participant_id %in% l_missing_fromV7)

setdiff(colnames(clin_data), colnames(clin_data_V7))
#[1] "registered_at_gmc_trust"              "clinic_sample_collected_at_gmc"       "clinic_sample_collected_at_gmc_trust"
#[4] "case_solved_family"                   "postdist"                             "rescty"                              
#[7] "superpopu"
clin_data_V7$registered_at_gmc_trust = rep("RE_V7",length(clin_data_V7$participant_id))
clin_data_V7$clinic_sample_collected_at_gmc = rep("RE_V7",length(clin_data_V7$participant_id))
clin_data_V7$clinic_sample_collected_at_gmc_trust = rep("RE_V7",length(clin_data_V7$participant_id))
clin_data_V7$case_solved_family = rep("RE_V7",length(clin_data_V7$participant_id))
clin_data_V7$postdist = rep("RE_V7",length(clin_data_V7$participant_id))
clin_data_V7$rescty = rep("RE_V7",length(clin_data_V7$participant_id))
clin_data_V7$superpopu = rep("RE_V7",length(clin_data_V7$participant_id))

# Reorder column names as `clin_data`
clin_data_V7 = clin_data_V7[names(clin_data)]

# Merge RE_V7 with clin_data
clin_data = rbind(clin_data,
                  clin_data_V7)
clin_data = unique(clin_data)
dim(clin_data)
#

# Let´s put all panel names into 1 single string splitted by ','
list_panels = clin_data %>% group_by(participant_id) %>% summarise(panel_list = toString(unique(panel_name))) %>% ungroup() %>% as.data.frame()
dim(list_panels)
# 87139  2

# Let´s put all HPO terms into 1 single string splitted by ','
list_hpos = clin_data %>% group_by(participant_id) %>% summarise(hpo_list = toString(unique(hpo_term))) %>% ungroup() %>% as.data.frame()
dim(list_hpos)
# 87139  2

# Let's put all specific diseases into 1 single string splitted by ','
list_diseases = clin_data %>% group_by(participant_id) %>% summarise(diseases_list = toString(unique(normalised_specific_disease))) %>% ungroup() %>% as.data.frame()
dim(list_diseases)
# 87139  2

list_disease_group = clin_data %>% group_by(participant_id) %>% summarise(diseasegroup_list = toString(unique(disease_group))) %>% ungroup() %>% as.data.frame()
dim(list_disease_group)
# 87139  2

list_disease_subgroup = clin_data %>% group_by(participant_id) %>% summarise(diseasesubgroup_list = toString(unique(disease_sub_group))) %>% ungroup() %>% as.data.frame()
dim(list_disease_subgroup)
# 87139  2

# Remove the panels and hpo columns, and include the list of panels and hpo respectively
clin_data = clin_data %>% 
  select(participant_id, platekey, rare_diseases_family_id, biological_relationship_to_proband, genetic_vs_reported_results, participant_ethnic_category, genome_build, year_of_birth, 
         participant_phenotypic_sex, participant_type, programme, type, family_group_type, affection_status, superpopu, clinic_sample_collected_at_gmc, clinic_sample_collected_at_gmc_trust, case_solved_family,
         registered_at_gmc_trust, rescty, postdist)
dim(clin_data)
# 2440099  21

clin_data = left_join(clin_data,
                      list_diseases,
                      by = "participant_id")
dim(clin_data)
# 2440099  22

clin_data = left_join(clin_data,
                      list_disease_group,
                      by = "participant_id")
dim(clin_data)
# 2440099 23

clin_data = left_join(clin_data,
                      list_disease_subgroup,
                      by = "participant_id")
dim(clin_data)
# 2440099 24

clin_data = left_join(clin_data,
                      list_panels,
                      by = "participant_id")
dim(clin_data)
# 2440099 25

clin_data = left_join(clin_data,
                      list_hpos,
                      by = "participant_id")
dim(clin_data)
# 2440099 26

# Enrich clin_data with pilot_clin_data, keeping diff fields as `.`
colnames(pilot_clin_data) = c("participant_id", "platekey", "rare_diseases_family_id", "participant_phenotypic_sex", "biological_relationship_to_proband", "affection_status", "year_of_birth", 
                              "ageOfOnset", "genetic_vs_reported_results", "diseases_list", "best_guess_predicted_ancstry", "bestGUESS_super_pop", "participant_ethnic_category", 
                              "diseasesubgroup_list", "diseasegroup_list")

# Generate extra columns from clin data for pilot clin data
pilot_clin_data$genome_build = rep("GRCh37", length(pilot_clin_data$participant_id))
pilot_clin_data$programme = rep("RD Pilot", length(pilot_clin_data$participant_id))
pilot_clin_data$family_group_type = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$panel_list = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$hpo_list = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$clinic_sample_collected_at_gmc = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$clinic_sample_collected_at_gmc_trust = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$case_solved_family = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$registered_at_gmc_trust = rep("Pilot", length(pilot_clin_data$participant_id))
pilot_clin_data$type = rep("rare disease germline", length(pilot_clin_data$participant_id))
pilot_clin_data$rescty = rep("Pilot", length(pilot_clin_data$participant_id))
pilot_clin_data$postdist = rep("Pilot", length(pilot_clin_data$participant_id))
# generate `participant_type` 
pilot_clin_data = pilot_clin_data %>%
  group_by(participant_id) %>%
  mutate(participant_type = ifelse(biological_relationship_to_proband %in% "Proband", "Proband", "Relative")) %>%
  ungroup() %>%
  unique() %>%
  as.data.frame()
  
# Select columns
pilot_clin_data = pilot_clin_data %>%
  select(participant_id, platekey, rare_diseases_family_id, biological_relationship_to_proband,
         genetic_vs_reported_results, participant_ethnic_category ,genome_build, year_of_birth, 
         participant_phenotypic_sex, participant_type, programme, type, family_group_type, affection_status, 
         bestGUESS_super_pop, clinic_sample_collected_at_gmc, clinic_sample_collected_at_gmc_trust, case_solved_family,
         registered_at_gmc_trust, rescty, postdist,
         diseases_list, diseasegroup_list, diseasesubgroup_list, panel_list, hpo_list)
colnames(pilot_clin_data) = colnames(clin_data)

clin_data = rbind(clin_data,
                  pilot_clin_data)
dim(clin_data)
# 2444933  26

# Check genQA cases
df_genQA = read.csv("~/Documents/STRs/VALIDATION/genQA/genQA/merged_batch1_batch3_STRs.tsv",
                    stringsAsFactors = F, 
                    header = T,
                    sep = "\t")
dim(df_genQA)
# 51  17

# Include genQA data in clin_data
df_genQA$programme = rep("genQA", length(df_genQA$Platekey))
df_genQA$type = rep("rare disease germline", length(df_genQA$Platekey))
df_genQA$biological_relationship_to_proband = rep("genQA", length(df_genQA$Platekey))
df_genQA$participant_type = rep("genQA", length(df_genQA$Platekey))
df_genQA$genetic_vs_reported_results = rep("genQA", length(df_genQA$Platekey))
df_genQA$participant_ethnic_category = rep("genQA", length(df_genQA$Platekey))
df_genQA$genome_build = rep("genQA", length(df_genQA$Platekey))
df_genQA$year_of_birth = rep("genQA", length(df_genQA$Platekey))
df_genQA$participant_phenotypic_sex = rep("genQA", length(df_genQA$Platekey))
df_genQA$family_group_type = rep("genQA", length(df_genQA$Platekey))
df_genQA$affection_status = rep("genQA", length(df_genQA$Platekey))
df_genQA$superpopu = rep("genQA", length(df_genQA$Platekey))
df_genQA$clinic_sample_collected_at_gmc = rep("genQA", length(df_genQA$Platekey))
df_genQA$clinic_sample_collected_at_gmc_trust = rep("genQA", length(df_genQA$Platekey))
df_genQA$case_solved_family = rep("genQA", length(df_genQA$Platekey))
df_genQA$diseases_list = rep("genQA", length(df_genQA$Platekey))
df_genQA$diseasegroup_list = rep("genQA", length(df_genQA$Platekey))
df_genQA$diseasesubgroup_list = rep("genQA", length(df_genQA$Platekey))
df_genQA$panel_list = rep("genQA", length(df_genQA$Platekey))
df_genQA$hpo_list = rep("genQA", length(df_genQA$Platekey))
df_genQA$registered_at_gmc_trust = rep("genQA", length(df_genQA$Platekey))
df_genQA$rescty = rep("genQA", length(df_genQA$Platekey))
df_genQA$postdist = rep("genQA", length(df_genQA$Platekey))

df_genQA = df_genQA %>%
  select(PID, Platekey, familyID, biological_relationship_to_proband,
         genetic_vs_reported_results, participant_ethnic_category, genome_build, year_of_birth,
         participant_phenotypic_sex, participant_type, programme, type, family_group_type, affection_status,
         superpopu, clinic_sample_collected_at_gmc, clinic_sample_collected_at_gmc_trust, case_solved_family,
         registered_at_gmc_trust, rescty, postdist,
         diseases_list, diseasegroup_list, diseasesubgroup_list, panel_list,
         hpo_list)
colnames(df_genQA) = colnames(clin_data)

clin_data = rbind(clin_data,
                  df_genQA)
dim(clin_data)
# 2444984  26

# Write into a file
write.table(clin_data,
            "Main_RE_V11_and_Pilot_programmes.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
