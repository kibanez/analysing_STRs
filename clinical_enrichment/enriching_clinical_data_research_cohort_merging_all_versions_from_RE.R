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


clin_data_merged = rbind(clin_data_merged,
                         clin_data_v4 %>% filter(participant_id %in% setdiff(pids_clin_data_v4, clin_data_merged$participant_id)))


l_pids = c(l_pids,
           setdiff(pids_clin_data_v4, l_pids))
length(l_pids)
# 1182293


# Merge with V3
pids_clin_data_v3 = read.csv("../clinical_data/raw/RE_clinical_data_V3/genome_file_paths_and_types_2020-07-07_11-17-23.tsv",
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
pids_clin_data_v3 = pids_clin_data_v3$participant_id

# No extra genomes
l_pids = c(l_pids,
           setdiff(pids_clin_data_v3, l_pids))
length(l_pids)
# 1182293

# Merge with V2
pids_clin_data_v2 = read.csv("../clinical_data/raw/RE_clinical_data_V2/participant_2020-07-07_11-23-34.tsv",
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
pids_clin_data_v2 = pids_clin_data_v2$participant_id

# 3,211 extra genomes
l_pids = c(l_pids,
           setdiff(pids_clin_data_v2, l_pids))
length(l_pids)
# 1185504

# Merge with V1
pids_clin_data_v1 = read.csv("../clinical_data/raw/RE_clinical_data_V1/genome_2020-07-07_11-27-22.tsv",
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
pids_clin_data_v1 = pids_clin_data_v1$participant_id

# 1 extra genome
l_pids = c(l_pids,
           setdiff(pids_clin_data_v1, l_pids))
length(l_pids)
# 1185505


clin_data_merged = clin_data_v9 %>%
  filter(participant_id %in% l_pids)
dim(clin_data_merged)
# 1180803  31

clin_data_merged = rbind(clin_data_merged,
                         clin_data_v8 %>% filter(participant_id %in% l_pids, !l_pids %in% clin_data_v9$participant_id))
  
  
clin_data_merged = unique(clin_data_merged)
dim(clin_data_merged)
# 

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
