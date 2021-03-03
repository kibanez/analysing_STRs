# Objective: from Matteo/Ari work analysing expansions on AR (from EHv255)
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/AR/feb2021/gvcfgenotyper/")

# Load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V10_and_Pilot_programmes.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 2101385  24

# Load unrel genomes (batch2)
l_unrel = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/l_unrelated_55603_genomes_batch2.txt",
                     stringsAsFactors = F)
l_unrel = l_unrel$V1
length(l_unrel)
# 55603

# Load 38 expanded genomes in AR- after visual QC inspection
l_exp_genomes = read.table("./list_38_genomes_beyond_patho.txt", stringsAsFactors = F)
l_exp_genomes = l_exp_genomes$V1
length(l_exp_genomes)
# 38 

# Retrieve genome build, sex, popu for these genomes
haplo_genomes = clin_data %>%
  filter(platekey %in% l_exp_genomes) %>%
  select(platekey, genome_build, participant_ethnic_category, participant_phenotypic_sex, superpopu, programme)
haplo_genomes = unique(haplo_genomes)
dim(haplo_genomes)
# 41  6

# There is 1 genome not included in the clinical data table, which belongs to data release V9
setdiff(l_exp_genomes,unique(haplo_genomes$platekey))
# "LP3000458-DNA_D11"

#to_add = c("LP3000458-DNA_D11", "GRCh38", "White: British", "Male", "EUR", "Rare diseases germline", "/genomes/by_date/2018-01-10/HX10655850/LP3000458-DNA_D11")
to_add = c("LP3000458-DNA_D11", "GRCh38", "White: British", "Male", "EUR", "Rare diseases germline")
haplo_genomes = rbind(haplo_genomes,
                      to_add)
dim(haplo_genomes)
# 42  6

# Create 2 groups: Females and Males
haplo_genomes_female = haplo_genomes %>%
  filter(participant_phenotypic_sex %in% "Female")
  
haplo_genomes_male = haplo_genomes %>%
  filter(participant_phenotypic_sex %in% "Male")


length(unique(haplo_genomes_female$platekey))
# 18
length(unique(haplo_genomes_male$platekey))
# 20

write.table(haplo_genomes_female,
            "./table_18genomes_female.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

write.table(haplo_genomes_male,
            "./table_20genomes_male.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# selecting genomes for CONTROL cohort
#  genome_build %in% GRCh38 , affection_status %in% unaffected, repeat_size < 37, and population %in% EUR since the majority of genomes in the cases cohort correspond to Europeans. 
# Unrelated, belonging to different families each of them.
# reported vs genetic checks PASS
table_ehv2 = read.csv("~/Documents/STRs/ANALYSIS/haplotyping/AR/HaploView/Matteo_Ari/table_STR_repeat_size_each_row_allele_EHv2.5.5_AR_CAG_simplified_dedup_050220.tsv",
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")
dim(table_ehv2)
# 115772  19

list_expanded = table_ehv2 %>% 
  filter(repeat_size >=37) %>%
  select(platekey) %>%
  pull() %>%
  unique()
length(list_expanded)
# 90

table_ehv2 = left_join(table_ehv2,
                       clin_data %>% select(platekey, genetic_vs_reported_results),
                       by = "platekey")

# Only proband is enriched with genetic_vs_reported checks. Let's take the famIDs that have been undergone successfully
# ONly rare diseases. We don't do this with cancer
l_famID_passesGvsRChecks = table_ehv2 %>%
  filter(genetic_vs_reported_results %in% "familyPassesGvsRChecks") %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull()

female_controls = table_ehv2 %>%
  filter(genome_build %in% "GRCh38", affection_status %in% "Unaffected", population %in% "EUR", !platekey %in% list_expanded, 
         participant_phenotypic_sex %in% "Female", rare_diseases_family_id %in% l_famID_passesGvsRChecks)
female_controls = unique(female_controls)
dim(female_controls)
# 17929  20

male_controls = table_ehv2 %>%
  filter(genome_build %in% "GRCh38", affection_status %in% "Unaffected", population %in% "EUR", !platekey %in% list_expanded,
         participant_phenotypic_sex %in% "Male", rare_diseases_family_id %in% l_famID_passesGvsRChecks)
male_controls = unique(male_controls)
dim(male_controls)
# 7308  20

# Filter out related genomes
female_controls = female_controls %>%
  filter(platekey %in% l_unrel)
dim(female_controls)
# 14596  20

male_controls = male_controls %>%
  filter(platekey %in% l_unrel)
dim(male_controls)
# 6705  20 

l_controls_AR_females = unique(female_controls$platekey)
length(l_controls_AR_females)
# 7709
l_controls_AR_males = unique(male_controls$platekey)
length(l_controls_AR_males)
# 6627

# Write into files whole list of male and female controls for AR
write.table(female_controls,
            "./table_female_7709_genomes_CONTROL_for_AR.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")
write.table(male_controls,
            "./table_male_6627_genomes_CONTROL_for_AR.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")

# Take path to the genome VCF files and write them into a file
upload_report = read.csv("~/Documents/STRs/clinical_data/clinical_data/upload_report.2020-08-18.txt",
                         stringsAsFactors = F,
                         sep = "\t",
                         header = T)
dim(upload_report)
# 120711  10

# we need to remove `_copied` - IT people changed again
upload_report$Platekey = gsub("_copied", "", upload_report$Platekey)
upload_report$Path = gsub("_copied", "", upload_report$Path)

# create new path
upload_report = upload_report %>%
  group_by(Platekey) %>%
  mutate(gvcf_path = paste(paste(paste(paste(Path, "Variations", sep = "/")),Platekey,sep = "/"),".genome.vcf.gz", sep = "")) %>%
  ungroup() %>%
  as.data.frame()

# Some platekeys have been sequenced more than once, let's select the latest one
upload_report = upload_report %>%
  group_by(Platekey) %>%
  mutate(latest_gvcf_path = max(gvcf_path)) %>%
  ungroup() %>%
  as.data.frame()

# Retrieve gVCF files for male AR cases that have been sequenced in GRCh38
list_gvcf_cases_male = upload_report %>%
  filter(Platekey %in% unique(haplo_genomes_male$platekey), Delivery.Version %in% "V4") %>%
  select(latest_gvcf_path) %>%
  unique() %>%
  pull()
length(list_gvcf_cases_male)
# 17

write.table(list_gvcf_cases_male,
            "list_17_gVCF_AR_male_CASES_GRCh38.txt",
            quote = F, row.names = F, col.names = F)

# Retrieve gVCF files for male AR controls that have been sequenced in GRCh38
list_gvcf_male = upload_report %>%
  filter(Platekey %in% l_controls_AR_males, Delivery.Version %in% "V4") %>%
  select(latest_gvcf_path) %>%
  unique() %>%
  pull()
length(list_gvcf_male)
# 6627

write.table(list_gvcf_male,
            "list_6627_gVCF_AR_male_CONTROLS_GRCh38.txt",
            quote = F, row.names = F, col.names = F)

# Retrieve gVCF files for female AR controls that have been sequenced in GRCh38
list_gvcf_female = upload_report %>%
  filter(Platekey %in% l_controls_AR_females, Delivery.Version %in% "V4") %>%
  select(latest_gvcf_path) %>%
  unique() %>%
  pull()
length(list_gvcf_female)
# 7709

write.table(list_gvcf_female,
            "list_7709_gVCF_AR_female_CONTROLS_GRCh38.txt",
            quote = F, row.names = F, col.names = F)

# Create phenotype file for males --> for plink
# FID IID PID MID sex affection
list_total_samples = read.table("cases_and_controls/list_16_cases_6627_controls_males_AR.txt", stringsAsFactors = F)
list_total_samples = list_total_samples$V1
length(list_total_samples)
# 6643

df_phenotype = clin_data %>%
  filter(platekey %in% list_total_samples) %>%
  select(platekey, participant_phenotypic_sex) %>%
  unique()
length(unique(df_phenotype$platekey))
# 6642

# Enrich with affection status following plink format
#-9 missing 
#0 missing
#1 unaffected
#2 affected
df_phenotype = df_phenotype %>%
  group_by(platekey) %>%
  mutate(affection = ifelse(platekey %in% l_exp_genomes, "2", "1")) %>%
  ungroup() %>%
  as.data.frame()

# "LP3000458-DNA_D11" is missing, male and expanded
df_phenotype = rbind(df_phenotype,
                     c("LP3000458-DNA_D11", "Male", "2"))
length(unique(df_phenotype$platekey))
# 6643


# Reformat gender (male to 1, female to 2)
df_phenotype$participant_phenotypic_sex = gsub("Male", "1", df_phenotype$participant_phenotypic_sex)
df_phenotype$participant_phenotypic_sex = gsub("Female", "2", df_phenotype$participant_phenotypic_sex)

# Include FID, PID, MID
df_phenotype$FID = df_phenotype$platekey
df_phenotype$PID = rep(0, length(df_phenotype$platekey))
df_phenotype$MID = rep(0, length(df_phenotype$platekey))

# Sort phenotype file
df_phenotype = df_phenotype %>%
  select(FID, platekey, PID, MID, participant_phenotypic_sex, affection)
colnames(df_phenotype) = c("FID", "IID", "PID", "MID", "sex", "affection")
write.table(df_phenotype, "cases_and_controls/phenotype_file.tsv", quote = F, col.names = T, row.names = F, sep = "\t")
