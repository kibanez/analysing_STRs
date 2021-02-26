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

l_controls_AR_females = unique(female_controls$platekey)
length(l_controls_AR_females)
# 9442
l_controls_AR_males = unique(male_controls$platekey)
length(l_controls_AR_males)
# 7219

# Write into files whole list of male and female controls for AR
write.table(female_controls,
            "./table_female_9442_genomes_CONTROL_for_AR.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")
write.table(male_controls,
            "./table_male_7219_genomes_CONTROL_for_AR.tsv",
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

# Retrieve gVCF files for male AR controls that have been sequenced in GRCh38
list_gvcf_male = upload_report %>%
  filter(Platekey %in% l_controls_AR_males, Delivery.Version %in% "V4") %>%
  select(latest_gvcf_path) %>%
  unique() %>%
  pull()
length(list_gvcf_male)
# 10617

write.table(list_gvcf_male,
            "list_10617_gVCF_AR_male_CONTROLS_GRCh38.txt",
            quote = F, row.names = F, col.names = F)

# Retrieve gVCF files for female AR controls that have been sequenced in GRCh38
list_gvcf_female = upload_report %>%
  filter(Platekey %in% l_controls_AR_females, Delivery.Version %in% "V4") %>%
  select(latest_gvcf_path) %>%
  unique() %>%
  pull()
length(list_gvcf_female)
# 13807

write.table(list_gvcf_female,
            "list_13807_gVCF_AR_female_CONTROLS_GRCh38.txt",
            quote = F, row.names = F, col.names = F)


# defining a specific seed, so every time we run this script, we end up selecting the "same random" genomes
set.seed(5)

# Take random 20 genomes per each male and female control groups
# Taking unrelated ones (by FamilyId)
l_random_20_family_female = sample(female_controls$rare_diseases_family_id, 20)
l_random_20_family_male = sample(male_controls$rare_diseases_family_id, 20)

random_20_female = female_controls %>%
  filter(rare_diseases_family_id %in% l_random_20_family_female) %>%
  select(rare_diseases_family_id, platekey, genome_build, participant_phenotypic_sex, population)

random_20_male = male_controls %>%
  filter(rare_diseases_family_id %in% l_random_20_family_male) %>%
  select(rare_diseases_family_id, platekey, genome_build, participant_phenotypic_sex, population)

# Enrich with the path to the gVCF
upload_report = read.csv("~/Documents/STRs/data/research/input/batch_august2020_EHv255_and_EHv322/input/upload_report.2020-08-18.txt",
                         stringsAsFactors = F,
                         sep = "\t",
                         comment.char = "#",
                         header = F)
# Remove `_copied` from V3
upload_report$V3 = gsub("_copied", "", upload_report$V3)

random_20_female = left_join(random_20_female,
                             upload_report %>% select(V3,V6),
                             by = c("platekey" = "V3"))

random_20_female$V6 = gsub("_copied", "/Variations/", random_20_female$V6)
random_20_female = random_20_female %>%
  group_by(platekey) %>%
  mutate(V6 = paste(V6, paste(platekey, ".genome.vcf.gz", sep = ""), sep = "")) %>%
  ungroup() %>%
  as.data.frame()
random_20_female = unique(random_20_female)
dim(random_20_female)
# 23 6

random_20_male = left_join(random_20_male,
                           upload_report %>% select(V3,V6),
                           by = c("platekey" = "V3"))
random_20_male$V6 = gsub("_copied", "/Variations/", random_20_male$V6)
random_20_male = random_20_male %>%
  group_by(platekey) %>%
  mutate(V6 = paste(V6, paste(platekey, ".genome.vcf.gz", sep = ""), sep = "")) %>%
  ungroup() %>%
  as.data.frame()
random_20_male = unique(random_20_male)
dim(random_20_male)
# 20 6

# Write them into files
write.table(random_20_female,
            "./table_20genomes_CONTROL_female.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")

write.table(random_20_male,
            "./table_20genomes_CONTROL_male.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")
