# The objective here is to enrich with clinical data all platekeys/participantIDs we do have selected for the research cohort study
# Merging all versions RE has published so far (not to lose anything)
# Joining Main and Pilot 
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr)

# setwd
setwd("~/Documents/STRs/clinical_data/clinical_research_cohort/")

# Clinical data retrieved from RE, diff versions
# Complete until V4
# V1,V2,V3 - take the list of genomes, and check/confirm they are all included in any of the versions from V4

# NOT RUN ANYMORE -- go to line 97
clin_data_v10 = read.csv("../clinical_data/rd_genomes_all_data_080920_V10.tsv",
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = "\t")
dim(clin_data_v10)
# 1184730  31

clin_data_v9 = read.csv("../clinical_data/rd_genomes_all_data_100720_V9.tsv",
                         header = TRUE,
                         stringsAsFactors = FALSE,
                         sep = "\t")
dim(clin_data_v9)
# 1180803  31

clin_data_v8 = read.csv("../clinical_data/rd_genomes_all_data_100720_V8.tsv",
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = "\t")
dim(clin_data_v8)
# 1124633  31

clin_data_v7 = read.csv("../clinical_data/rd_genomes_all_data_100720_V7.tsv",
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = "\t")
dim(clin_data_v7)
# 1021801  31

clin_data_v6 = read.csv("../clinical_data/rd_genomes_all_data_100720_V6.tsv",
                       header = TRUE,
                       stringsAsFactors = FALSE,
                       sep = "\t")
dim(clin_data_v6)
# 775234  31

clin_data_v5 = read.csv("../clinical_data/rd_genomes_all_data_100720_V5.1.tsv",
                       header = TRUE,
                       stringsAsFactors = FALSE,
                       sep = "\t")
dim(clin_data_v5)
# 555656  31

# Let's merge V5-V10 version tables, but keeping the info frmo latest version if PID is the same (there might be inconsistencies)
# Check and confirm V4 has not extra genomes
# Check and confirm V1-V3 have not extra genomes
colnames(clin_data_v10) = colnames(clin_data_v9)
colnames(clin_data_v8) = colnames(clin_data_v9)
colnames(clin_data_v7) = colnames(clin_data_v9)
colnames(clin_data_v6) = colnames(clin_data_v9)
colnames(clin_data_v5) = colnames(clin_data_v9)

# V10 with V9
clin_data_merged = clin_data_v10

pids_clin_data_v9 = unique(clin_data_v9$participant_id)
l_new_pid_v9 = setdiff(pids_clin_data_v9, clin_data_merged$participant_id)
length(l_new_pid_v9)
# 272

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v9 %>% filter(participant_id %in% l_new_pid_v9))
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 1168468  31

# V9 with V8
pids_clin_data_v8 = unique(clin_data_v8$participant_id)
l_new_pid_v8 = setdiff(pids_clin_data_v8, clin_data_merged$participant_id)
length(l_new_pid_v8)
# 331

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v8 %>% filter(participant_id %in% l_new_pid_v8))
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 1172993 31

# merged with V7
pids_clin_data_v7 = unique(clin_data_v7$participant_id)
l_new_pid_v7 = setdiff(pids_clin_data_v7, clin_data_merged$participant_id)
length(l_new_pid_v7)
# 865

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v7 %>% filter(participant_id %in% l_new_pid_v7))
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 1186715 31

# merged with V6
pids_clin_data_v6 = unique(clin_data_v6$participant_id)
l_new_pid_v6 = setdiff(pids_clin_data_v6, clin_data_merged$participant_id)
length(l_new_pid_v6)
# 271 

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v6 %>% filter(participant_id %in% l_new_pid_v6))
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 1189208 31

# merged with V5
pids_clin_data_v5 = unique(clin_data_v5$participant_id)
l_new_pid_v5 = setdiff(pids_clin_data_v5, clin_data_merged$participant_id)
length(l_new_pid_v5)
# 10 

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v5 %>% filter(participant_id %in% l_new_pid_v5))
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 1189268 31

# Let's write this into file, not to re-run
write.table(clin_data_merged,
            "./clin_data_merged_V5:V10.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")

# Load from the beginning all merged data
clin_data_merged = read.csv("./clin_data_merged_V5:V10.tsv",
                            stringsAsFactors = F,
                            header = T,
                            sep = "\t")
dim(clin_data_merged)
# 1189268 31

# Merge with V4 
clin_data_v4 = read.csv("../clinical_data/rd_genomes_all_data_100720_V4.tsv",
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = "\t")
dim(clin_data_v4)
# 370296  29

pids_clin_data_v4 = read.csv("../clinical_data/raw/RE_clinical_data_V4/genome_file_paths_and_types_2020-07-07_10-55-46.tsv",
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
pids_clin_data_v4 = pids_clin_data_v4$participant_id

# there is 1 extra genome
length(setdiff(pids_clin_data_v4, clin_data_merged$participant_id))
# 1
setdiff(pids_clin_data_v4, clin_data_merged$participant_id)
# 115002748

# clin_data_v4 is missing 2 columns clin_data_merged has: `affection_status` and `genetic_vs_reported_results`
clin_data_v4 = clin_data_v4 %>%
  filter(participant_id %in% setdiff(pids_clin_data_v4, clin_data_merged$participant_id)) %>%
  unique()
clin_data_v4$affection_status = '.'
clin_data_v4$genetic_vs_reported_results = '.'

clin_data_v4 = clin_data_v4 %>%
  select(participant_id, platekey.x, file_path, type, rare_diseases_family_id, platekey.y, biological_relationship_to_proband, participant_type,
         normalised_specific_disease, genome_build, genetic_vs_reported_results, participant_ethnic_category, disease_group, disease_sub_group,
         specific_disease, participant_medical_review_qc_state_code, year_of_birth, participant_phenotypic_sex, participant_karyotypic_sex, participant_stated_gender,
         programme_consent_status, programme, family_group_type, family_medical_review_qc_state_code, panel_name, panel_version, hpo_term, hpo_id, affection_status, 
         best_guess_predicted_ancstry, self_reported)
colnames(clin_data_v4) = colnames(clin_data_merged)

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v4)
dim(clin_data_merged)
# 1180678  31

# Merge with V3
pids_clin_data_v3 = read.csv("../clinical_data/raw/RE_clinical_data_V3/genome_file_paths_and_types_2020-07-07_11-17-23.tsv",
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
pids_clin_data_v3 = pids_clin_data_v3$participant_id

# No extra genomes
length(setdiff(pids_clin_data_v3,
               clin_data_merged$participant_id))
# 0

# Merge with V2
pids_clin_data_v2 = read.csv("../clinical_data/raw/RE_clinical_data_V2/participant_2020-07-07_11-23-34.tsv",
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
pids_clin_data_v2 = unique(pids_clin_data_v2$participant_id)
length(pids_clin_data_v2)
# 53190

# 3,211 extra genomes
length(setdiff(pids_clin_data_v2,
               clin_data_merged$participant_id))
# 3211

clin_data_v2 = read.csv("../clinical_data/rd_genomes_all_data_130720_V2.tsv",
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = "\t")
dim(clin_data_v2)
# 100933  26

#include not existing fields
clin_data_v2$genetic_vs_reported_results = rep('.', length(clin_data_v2$participant_id))
clin_data_v2$panel_name = rep('.', length(clin_data_v2$participant_id))
clin_data_v2$panel_version = rep('.', length(clin_data_v2$participant_id))
clin_data_v2$affection_status = rep('.', length(clin_data_v2$participant_id))
clin_data_v2$plate_key = clin_data_v2$platekey
clin_data_v2 = clin_data_v2 %>%
  select(participant_id, platekey, path, type, rare_diseases_family_id, plate_key, biological_relationship_to_proband, participant_type,
         normalised_specific_disease, genome_build, genetic_vs_reported_results, participant_ethnic_category,
         disease_group, disease_sub_group, specific_disease, participant_medical_review_qc_state_code,
         year_of_birth, participant_phenotypic_sex, participant_karyotypic_sex, participant_stated_gender,
         programme_consent_status, programme, family_group_type, family_medical_review_qc_state_code,
         panel_name, panel_version, hpo_term, hpo_id, 
         affection_status, best_guess_predicted_ancstry, self_reported)
colnames(clin_data_v2) = colnames(clin_data_merged)

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v2)
dim(clin_data_merged)
# 1281611  31

# Merge with V1
pids_clin_data_v1 = read.csv("../clinical_data/raw/RE_clinical_data_V1/genome_2020-07-07_11-27-22.tsv",
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
pids_clin_data_v1 = unique(pids_clin_data_v1$participant_id)
length(pids_clin_data_v1)
# 18446

# 1 extra genome
length(setdiff(pids_clin_data_v1,
               clin_data_merged$participant_id))
# 1
setdiff(pids_clin_data_v1,
        clin_data_merged$participant_id)
# 111001387

# let's take the info corresponding to this PID from the RE V1 raw data
line_v1 = c(setdiff(pids_clin_data_v1,
                    clin_data_merged$participant_id),
            "LP3000143-DNA_F04",
            "/genomes/by_date/2016-12-06/HX01615322/LP3000143-DNA_F04/Assembly/LP3000143-DNA_F04.bam",
            "rare disease germline",
            "111001380",
            "LP3000143-DNA_F04",
            "Mother",
            "Relative",
            NA,
            "GRCh37",
            NA,
            "Not Stated",
            NA,
            NA,
            NA,
            "Passed medical review - for interpretation",
            NA,
            "Female",
            "Not Supplied",
            "Female",
            "Consenting",
            "Rare Diseases",
            "Trio with Mother and Father",
            "Passed medical review - for interpretation",
            NA,
            NA,
            NA,
            NA,
            "Unaffected",
            NA,
            NA)

clin_data_merged = rbind(clin_data_merged,
                         line_v1)
dim(clin_data_merged)
# 1281612  31

# Load withdrawn pids and tag them
l_withdrawn = read.table("../clinical_data/raw/l_withdrawn_pid.txt", stringsAsFactors = F)
l_withdrawn = l_withdrawn$V1  
length(l_withdrawn)
# 27

index_withdrawn = which(clin_data_merged$participant_id %in% l_withdrawn)
clin_data_merged$withdrawn = "No"
clin_data_merged$withdrawn[index_withdrawn] = "Yes"

# Include PILOT data
pilot_data = read.csv("../pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_withPanels_280919.tsv",
                      stringsAsFactors = F,
                      sep = "\t",
                      header = T)
dim(pilot_data)
# 4974  11

# Include file path
pilot_path = read.csv("../pilot_clinical_data/b37_genomes_2019-12-11_11-07-20.tsv",
                      stringsAsFactors = F,
                      sep = "\t",
                      header = T)
dim(pilot_path)
# 4828  5

pilot_data = left_join(pilot_data,
                       pilot_path %>% select(participant_id, plate_key, path),
                       by = c("gelID" = "participant_id"))

pilot_data$type = rep("Pilot", length(pilot_data$gelID))
pilot_data$plate_key = pilot_data$plateKey
pilot_data = pilot_data %>%
  group_by(gelID) %>%
  mutate(participant_type = case_when(
    biological_relation_to_proband == "Proband" ~ "Proband",
    biological_relation_to_proband != "Proband" ~ "Relative"
  )) %>%
  ungroup() %>%
  as.data.frame()
pilot_data$normalised_specific_disease = pilot_data$specificDisease
pilot_data$genome_build = rep("GRCh37", length(pilot_data$gelID))
pilot_data$participant_ethnic_category = rep(NA, length(pilot_data$gelID))

# Let's enrich `disease_group` and `disease_subgroup` with pilot clin data
pheno_table = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/phenotyping_v140_2019-09-13_15-26-02.tsv",
                       stringsAsFactors = F,
                       header = T,
                       sep = "\t")
pheno_table = unique(pheno_table)
dim(pheno_table)
# 106  3

pilot_data = left_join(pilot_data,
                       pheno_table,
                       by = c("panel_list" = "specific_disease"))
pilot_data = unique(pilot_data)
dim(pilot_data)
# 4974  20

# Enrich with popu data
popu_pilot = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                      stringsAsFactors = F)
dim(popu_pilot)
# 4821  44

pilot_data = left_join(pilot_data,
                       popu_pilot %>% select(ID, bestGUESS_sub_pop, bestGUESS_super_pop),
                       by = c("plateKey" = "ID"))
dim(pilot_data)
# 4976  22

pilot_data$participant_medical_review_qc_state_code = pilot_data$qc_state
pilot_data$participant_karyotypic_sex = rep(NA, length(pilot_data$gelID)) 
pilot_data$programme_consent_status = rep("Consenting", length(pilot_data$gelID))
pilot_data$programme = rep("Pilot", length(pilot_data$gelID))
pilot_data$family_group_type = rep(".", length(pilot_data$gelID))
pilot_data$family_medical_review_qc_state_code = pilot_data$qc_state
pilot_data$panel_version = rep(NA, length(pilot_data$gelID))
pilot_data$normalised_specific_disease = pilot_data$specificDisease

# Enrich with hpo terms
hpo_table = read.csv("../pilot_clinical_data/data_freeze_Pilot_LK_RESEARCH/hpo.csv",
                     stringsAsFactors = F)
dim(hpo_table)
# 62501  8

pilot_data$gelID = as.character(pilot_data$gelID)

pilot_data = left_join(pilot_data,
                       hpo_table %>% filter(termPresence %in% "true") %>% select(gelID, term),
                       by = "gelID")

pilot_data$hpo_id = rep(NA, length(pilot_data$gelID))
pilot_data$withdrawn = rep("No", length(pilot_data$gelID))
pilot_data$participant_stated_gender = pilot_data$sex

pilot_data = pilot_data %>%
  select(gelID, plateKey, path, type,
         gelFamilyId.x, plate_key, biological_relation_to_proband, participant_type,
         normalised_specific_disease, genome_build, qc_state, participant_ethnic_category,
         disease_group, disease_subgroup, specificDisease, participant_medical_review_qc_state_code,
         yearOfBirth, sex, participant_karyotypic_sex, participant_stated_gender,
         programme_consent_status, programme, family_group_type, family_medical_review_qc_state_code,
         panel_list, panel_version, hpo_id, term,
         disease_status, bestGUESS_sub_pop, bestGUESS_super_pop, withdrawn)

colnames(pilot_data) = colnames(clin_data_merged)

clin_data_merged = rbind(clin_data_merged,
                         pilot_data)
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 1299442  32

# Write this into a file
write.table(clin_data_merged,
            "./clin_data_merged_V1:V9.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")


# List of all VCF files we have generated by running EH, march 2020
l_vcf_EHv255 = read.table("~/Documents/STRs/data/research/batch_march2020/list_EHv255_vcfs.txt",
                          stringsAsFactors = F)
l_vcf_EHv255 = l_vcf_EHv255$V1
length(l_vcf_EHv255)
# 92665

l_vcf_EHv322 = read.table("~/Documents/STRs/data/research/batch_march2020/list_EHv322_vcfs.txt",
                          stringsAsFactors = F)
l_vcf_EHv322 = l_vcf_EHv322$V1
length(l_vcf_EHv322)
# 92663

l_vcf = union(l_vcf_EHv255,
              l_vcf_EHv322)
length(l_vcf)
# 92665

l_vcf = gsub("./EH_", "", l_vcf)
l_vcf = gsub(".vcf", "", l_vcf)

length(intersect(l_vcf, clin_data_merged$platekey))
# 92235
length(setdiff(l_vcf, clin_data_merged$platekey))
# 430

# I'll keep like this, only containing clinical data from RE, not including that data in Catalog
# Let's prepare the data
# each row = 1 genome with specific disease, dis group, dis subgroup, HPO, sex, year of birth, participant type, affection status, population (estimated), and self-reported ancestry
list_norm_disease = clin_data_merged %>% group_by(participant_id) %>% summarise(list_norm_disease = toString(unique(normalised_specific_disease))) %>% ungroup() %>% as.data.frame()
list_panels = clin_data_merged %>% group_by(participant_id) %>% summarise(list_panels = toString(unique(panel_name))) %>% ungroup() %>% as.data.frame()
list_panels_version = clin_data_merged %>% group_by(participant_id) %>% summarise(list_panels_version = toString(unique(panel_version))) %>% ungroup() %>% as.data.frame()
list_hpos = clin_data_merged %>% group_by(participant_id) %>% summarise(list_hpos = toString(unique(hpo_term))) %>% ungroup() %>% as.data.frame()
list_hpos_id = clin_data_merged %>% group_by(participant_id) %>% summarise(list_hpos_id = toString(unique(hpo_id))) %>% ungroup() %>% as.data.frame()
list_dis_group = clin_data_merged %>% group_by(participant_id) %>% summarise(list_disease_group = toString(unique(disease_group))) %>% ungroup() %>% as.data.frame()
list_dis_subgroup = clin_data_merged %>% group_by(participant_id) %>% summarise(list_disease_subgroup = toString(unique(disease_sub_group))) %>% ungroup() %>% as.data.frame()

clin_data_merged = left_join(clin_data_merged,
                             list_norm_disease,
                             by = "participant_id")

clin_data_merged = left_join(clin_data_merged,
                             list_panels,
                             by = "participant_id")

clin_data_merged = left_join(clin_data_merged,
                             list_panels_version,
                             by = "participant_id")

clin_data_merged = left_join(clin_data_merged,
                             list_hpos,
                             by = "participant_id")

clin_data_merged = left_join(clin_data_merged,
                             list_hpos_id,
                             by = "participant_id")

clin_data_merged = left_join(clin_data_merged,
                             list_dis_group,
                             by = "participant_id")


clin_data_merged = left_join(clin_data_merged,
                             list_dis_subgroup,
                             by = "participant_id")

# Remove old columns
drops <- c("normalised_specific_disease","panel_name","panel_version", "hpo_term", "hpo_id", "disease_group", "disease_sub_group")
clin_data_merged = clin_data_merged[ , !(names(clin_data_merged) %in% drops)]
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 152193  32

# Remove duplicates
l_all_pids = unique(clin_data_merged$participant_id)
length(l_all_pids)
# 93614

# There are many duplicates where genetic_vs_reported have diff values - let's take the ones not NA
dedup_clin_data_merged = filter(clin_data_merged, (!is.na(participant_medical_review_qc_state_code) & programme %in% "Rare Diseases") | programme %in% c("Pilot", "Cancer"))
dim(dedup_clin_data_merged)
# 152187  32

length(unique(dedup_clin_data_merged$participant_id))
# 93612

# There are 2 PIDs which might be missing something - they are missing `Cancer` as programme
aver = clin_data_merged %>% filter(participant_id %in% setdiff(l_all_pids, unique(dedup_clin_data_merged$participant_id)))
aver$programme = rep("Cancer", length(aver$participant_id))

dedup_clin_data_merged = rbind(dedup_clin_data_merged,
                               aver)
dim(dedup_clin_data_merged)
# 152191  32

length(unique(dedup_clin_data_merged$participant_id))
# 93614

# Let's filter out types = `cancer tumour``, `cancer somatic`, and `experimental somatic`
dedup_clin_data_merged = dedup_clin_data_merged %>%
  filter(!type %in% c("cancer tumour", "cancer somatic", "experimental somatic"))
dim(dedup_clin_data_merged)
# 131769  32

length(unique(dedup_clin_data_merged$participant_id))
# 93491

# Let's filter out those not passing genetics vs re (familyPassesGvsRChecks)
dedup_clin_data_merged = dedup_clin_data_merged %>%
  filter((genetic_vs_reported_results %in% "familyPassesGvsRChecks" & programme %in% "Rare Diseases") | programme %in% c("Cancer", "Pilot"))
dim(dedup_clin_data_merged)
# 88795  32

length(unique(dedup_clin_data_merged$participant_id))
# 82506

# Select fields to keep 
#specific disease, dis group, dis subgroup, HPO, sex, year of birth, participant type, affection status, population (estimated), and self-reported ancestry
dedup_clin_data_merged = dedup_clin_data_merged %>%
  select(participant_id, rare_diseases_family_id, biological_relationship_to_proband, family_group_type, list_norm_disease, list_disease_group, list_disease_subgroup, list_hpos, list_hpos_id, participant_phenotypic_sex,
         year_of_birth, participant_type, affection_status, best_guess_predicted_ancstry, self_reported, participant_ethnic_category, programme,
         withdrawn, programme_consent_status, )
dedup_clin_data_merged = unique(dedup_clin_data_merged)
dim(dedup_clin_data_merged)
# 85237 19

length(unique(dedup_clin_data_merged$participant_id))
# 82506

l_dups = unique(dedup_clin_data_merged$participant_id[which(duplicated(dedup_clin_data_merged$participant_id))])
length(l_dups)
# 2676

# There are some PIDs with GRCh37 and GRCh38 genomes
dup_table = dedup_clin_data_merged %>% 
  filter(participant_id %in% l_dups)
dim(dup_table)
# 5407  19
length(unique(dup_table$participant_id))
# 2676

dedup_clin_data_merged = dedup_clin_data_merged %>%
  filter(!participant_id %in% l_dups)
dim(dedup_clin_data_merged)
# 79830  19

# Select those having less NAs
dup_table2 = data.frame()
l_dup_pid = unique(dup_table$participant_id)
for (i in 1:length(l_dup_pid)){
  aux = dup_table %>% filter(participant_id %in% l_dup_pid[i])
  
  # select the row with less NAs
  row_fewest_NAs = aux[which.max(rowSums(!is.na(aux))),]
  dup_table2 = rbind(dup_table2,
                     row_fewest_NAs)
}
dim(dup_table2)
# 2676  19
length(unique(dup_table2$participant_id))
# 2676

dedup_clin_data_merged = rbind(dedup_clin_data_merged,
                               dup_table2)
dedup_clin_data_merged = unique(dedup_clin_data_merged)
dim(dedup_clin_data_merged)
# 82506  19

length(unique(dedup_clin_data_merged$participant_id))
# 82506

write.table(dedup_clin_data_merged, 
            file = "clinical_data_research_cohort_82506_PIDs_merging_RE_V1toV9.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)
