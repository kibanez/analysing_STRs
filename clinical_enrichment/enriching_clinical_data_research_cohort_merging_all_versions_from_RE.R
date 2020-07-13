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

# Let's merge V5-V9 version tables
# Check and confirm V4 has not extra genomes
# Check and confirm V1-V3 have not extra genomes
colnames(clin_data_v8) = colnames(clin_data_v9)
colnames(clin_data_v7) = colnames(clin_data_v9)
colnames(clin_data_v6) = colnames(clin_data_v9)
colnames(clin_data_v5) = colnames(clin_data_v9)

clin_data_merged = rbind(clin_data_v9,
                         clin_data_v8)
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 1178166  31

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v7)
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 1282056  31

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v6)
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 1514861  31

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v5)
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 2061403  31

# Let's write this into file, not to re-run
write.table(clin_data_merged,
            "./clin_data_merged_V5:V9.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")

# Load from the beginning all merged data
clin_data_merged = read.csv("./clin_data_merged_V5:V9.tsv",
                            stringsAsFactors = F,
                            header = T,
                            sep = "\t")
dim(clin_data_merged)
# 2061403  31

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
# 2061409  31

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
# 2162342  31

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
# 2162343  31

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
# 2180173  32


# List of all VCF files we have generated by running EH
l_vcf = read.table("./list_vcf_research_cohort.tsv", stringsAsFactors = F)
l_vcf = l_vcf$V1
length(l_vcf)
# 86457


# Separate in GRCh37 and GRCh38
clinical_data_b37 = clinical_data %>% filter(genome_build %in% "GRCh37")
dim(clinical_data_b37)
# 136800  26

clinical_data_b38 = clinical_data %>% filter(genome_build %in% "GRCh38")
dim(clinical_data_b38)
# 877169  26


clinical_data_na = clinical_data %>% filter(is.na(genome_build))
dim(clinical_data_na)
# 35898 26

# Create a table with clinical data 
#LP_number Gel_ID	Family_ID	Relationship 	Affection status	Gender	YearOfBirth	Disease	Disease_subgroup	Disease_group

research_b37 = clinical_data_b37 %>% 
  select(plate_key.x, participant_id, rare_diseases_family_id, biological_relationship_to_proband, affection_status, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, disease_sub_group, disease_group, panel_name)
dim(research_b37)
# 136800  11

research_b37  = research_b37 %>% 
  filter(plate_key.x %in% l_vcf)
dim(research_b37)
# 113037  11


research_b38 = clinical_data_b38 %>% 
  select(plate_key.x, participant_id, rare_diseases_family_id, biological_relationship_to_proband, affection_status, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, disease_sub_group, disease_group, panel_name)
dim(research_b38)
# 877169  11

research_b38 = research_b38 %>%
  filter(plate_key.x %in% l_vcf)
dim(research_b38)
#  766399  11


research_na = clinical_data_na %>%
  select(plate_key.x, participant_id, rare_diseases_family_id, biological_relationship_to_proband, affection_status, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, disease_sub_group, disease_group, panel_name)
dim(research_na)
# 35898  11

# Check whether ALL VCF files in l_vcf have been included and enriched in both `research_b37` and `research_b38` tables
check_vcf_b37 = research_b37 %>% select(plate_key.x) %>% unique() %>% pull()
check_vcf_b38 = research_b38 %>% select(plate_key.x) %>% unique() %>% pull()
check_vcf_na = research_na %>% select(plate_key.x) %>% unique() %>% pull()

check_vcf = c(check_vcf_b37,
              check_vcf_b38,
              check_vcf_na)
length(check_vcf)
# 103438

# Which are not enriched?? Cancer samples?? <-- yes they are, we should enrich them as type (RD or Cancer)
setdiff(l_vcf, check_vcf)

length(setdiff(l_vcf, check_vcf))
# 380

## OUTPUT FILE
# When writing the output file that will facilitate us the statistical outcome when assessing the importance of our study, we will print the following:
# LP_number	Gel_ID	Family_ID	Relationship 	Affection status	Gender	YearOfBirth	Disease	Disease_subgroup	Disease_group
# Since there are many as NA's, because they are cancer, we will finally put everything on the same table

# First let's put all panel names as a list separated by ','

to_write = clinical_data %>%
  select(plate_key.x, participant_id, programme, genome_build, programme_consent_status, rare_diseases_family_id, biological_relationship_to_proband, affection_status, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, disease_sub_group, disease_group, panel_name, family_group_type, family_medical_review_qc_state_code)
dim(to_write)
# 1056568  16

# Change all `NA` in `biological_relationship_to_proband` to "proband"
l_proband = which(to_write$biological_relationship_to_proband == "N/A")
to_write$biological_relationship_to_proband[l_proband] = "Proband"

list_panels = to_write %>% group_by(participant_id) %>% summarise(panel_list = toString(panel_name)) %>% ungroup() %>% as.data.frame()
dim(list_panels)
# 89163  2

# Remove `panel_name` from to_write
to_write = to_write[,-14]
dim(to_write)
# 1056568  15

to_write = left_join(to_write,
                     list_panels,
                     by = "participant_id")

to_write = unique(to_write)

to_write %>% select(plate_key.x) %>% unique() %>% pull() %>% length()
# 113696

write.table(to_write, 
            file = "clinical_data_research_cohort_113696_genomes_withPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

# Since the whole table is super big, let's separate it
to_write_b38 = to_write %>%
  filter(genome_build %in% "GRCh38")
dim(to_write_b38)
# 69943  16

to_write_b37_na = to_write %>%
  filter(!genome_build %in% "GRCh38")
dim(to_write_b37_na)
# 47290  16


write.table(to_write_b38, 
            file = "clinical_data_research_cohort_113696_genomes_GRCh38_withPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

write.table(to_write_b37_na, 
            file = "clinical_data_research_cohort_113696_genomes_GRCh37_NA_withPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)


# Let's write the same info but removing the list of panels
write.table(to_write[,-16], 
            file = "clinical_data_research_cohort_113696_genomes_removingPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

write.table(to_write_b38[,-16], 
            file = "clinical_data_research_cohort_113696_genomes_GRCh38_removingPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

write.table(to_write_b37_na[,-16], 
            file = "clinical_data_research_cohort_113696_genomes_GRCh37_NA_removingPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)



## ONLY OUTPUT FILES considering 86,457 genomes we have analysed through EHv2.5.5 and EHv3.0.0
clinical_data_research_genomes = to_write %>% filter(plate_key.x %in% l_vcf)
dim(clinical_data_research_genomes)
# 89090  16

clinical_data_research_genomes %>% select(plate_key.x) %>% unique() %>% pull() %>% length()
# 86242

# Let's annotate the ones not included there with '.'
l_not_included = setdiff(l_vcf, clinical_data_research_genomes %>% select(plate_key.x) %>% unique() %>% pull())
length(l_not_included)
# 215

rep_char = rep(".", 15)

df_not_included = as.data.frame(t(c(l_not_included[1], rep_char)), stringsAsFactors=F)
colnames(df_not_included) = colnames(clinical_data_research_genomes)
for(i in 2:length(l_not_included)){
  df_not_included = rbind(df_not_included,
                          c(l_not_included[i], rep_char))
}


clinical_data_research_genomes = rbind(clinical_data_research_genomes,
                                       df_not_included)
dim(clinical_data_research_genomes)
# 89305  16

clinical_data_research_genomes %>% select(plate_key.x) %>% unique() %>% pull() %>% length()
# 86457


# Let's write this with and without panel names

write.table(clinical_data_research_genomes, 
            file = "clinical_data_research_cohort_86457_genomes_withPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

clinical_data_research_genomes %>% filter(genome_build %in% "GRCh38") %>% select(plate_key.x) %>% unique() %>% pull() %>% length()
# 58737 genomes
clinical_data_research_genomes %>% filter(genome_build %in% "GRCh38") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 58434 participants

write.table(clinical_data_research_genomes %>% filter(genome_build %in% "GRCh38"), 
            file = "clinical_data_research_cohort_86457_genomes_GRCh38_withPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

clinical_data_research_genomes %>% filter(!genome_build %in% "GRCh38") %>% select(plate_key.x) %>% unique() %>% pull() %>% length()
# 29097
clinical_data_research_genomes %>% filter(!genome_build %in% "GRCh38") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 19041

write.table(clinical_data_research_genomes %>% filter(!genome_build %in% "GRCh38"), 
            file = "clinical_data_research_cohort_86457_genomes_GRCh37_and_NA_withPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)


# The same without panels

write.table(clinical_data_research_genomes[,-16], 
            file = "clinical_data_research_cohort_86457_genomes_removingPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

clinical_data_research_genomes %>% filter(genome_build %in% "GRCh38") %>% select(plate_key.x) %>% unique() %>% pull() %>% length()
# 58737
write.table(clinical_data_research_genomes[,-16] %>% filter(genome_build %in% "GRCh38"), 
            file = "clinical_data_research_cohort_86457_genomes_GRCh38_removingPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

clinical_data_research_genomes %>% filter(!genome_build %in% "GRCh38") %>% select(plate_key.x) %>% unique() %>% pull() %>% length()
# 29097
write.table(clinical_data_research_genomes[,-16] %>% filter(!genome_build %in% "GRCh38"), 
            file = "clinical_data_research_cohort_86457_genomes_GRCh37_and_NA_removingPanels_250919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)
