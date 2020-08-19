# Objective: select genomes (better deduplicated) in order to run afterwards EH through all them
# The idea is to merge or join all data available to us: RE batch + population aggregated gVCF (in case we have more genomes to fish) + earlier EH batches
# RE batch (V9 latest batch from April 2020) - but we do have now the UNION of all RE batches (V1:V9)
# Population aggregated gVCF file (batch1 + batch2)
# EHv2 batch from march 2020
# EHv3 batch from march 2020

# Workflow: we will use Population aggregated gVCF file, as well as earlier EHv2 and EHv3 batches to complete the data
# But we will start from all genomes/platekeys sequenced so far at GEL, that have participantId info in Catalog
# Catalog studies: RD b38, RD b37, and Cancer (they have been selected by taking the info from the latest cohort)
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string # "R version 3.6.2 (2019-12-12)"

# Libraries
library(dplyr)

# Working directory
setwd("~/Documents/STRs/data/research/input")

# Load table with all genomes together with their paths
df_all_genomes = read.csv("~/Documents/STRs/data/research/input/batch_august2020_EHv255_and_EHv322/upload_report.2020-08-18.txt",
                          sep = "\t",
                          stringsAsFactors = F,
                          comment.char = "#",
                          header = F)
dim(df_all_genomes)
# 120711  2

# Load all merged clinical data from RE
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_research_cohort/clin_data_merged_V5:V9.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 2061403  31

# Let's keep only with the important information
clin_data = clin_data %>%
  select(platekey, participant_id, participant_phenotypic_sex, genome_build)
clin_data = unique(clin_data)
dim(clin_data)
# 170150 4

pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           sep = "\t",
                           stringsAsFactors = F,
                           header = T)
dim(pilot_clin_data)
# 4974  10

pilot_clin_data = pilot_clin_data %>%
  select(plateKey, gelID, sex)
dim(pilot_clin_data)
# 4974  3

pilot_clin_data$genome_build = rep("GRCh37", length(pilot_clin_data$plateKey))
colnames(pilot_clin_data) = colnames(clin_data)

# merge MAIN and PILOT clin data datasets
merged_clin_data = rbind(clin_data,
                         pilot_clin_data)
dim(merged_clin_data)
# 175124  4


# define column names
colnames(catalog_rd_b37) = c("cohort_id", "platekey", "participant_id", "isProband", "karyo", "sex", "affection_status", "build", "programme")
colnames(catalog_rd_b38) = c("cohort_id", "platekey", "participant_id", "isProband","sex", "affection_status", "build", "programme")

# remove those cohorts having REVOKED ord DEPRECATED
catalog_rd_b38 = catalog_rd_b38 %>% filter(!grepl("REVOKED", cohort_id))
catalog_rd_b38 = catalog_rd_b38 %>% filter(!grepl("DEPRECATED", cohort_id))

# Since we selected from Catalog info from THE LATEST COHORT DEFINED, we can ignore the cohort Id
catalog_rd_b37 = catalog_rd_b37 %>% 
  select(platekey, participant_id, build, programme)
catalog_rd_b37 = unique(catalog_rd_b37)
dim(catalog_rd_b37)
# 99  4

catalog_rd_b38 = catalog_rd_b38 %>%
  select(platekey, participant_id, build, programme)
catalog_rd_b38 = unique(catalog_rd_b38)
dim(catalog_rd_b38)
# 76881  4

# Cancer -> we will take this from the research environemnt

# Check duplicated info. Take deduplicated participantIDs together with their platekeys
unique_pid_b37 = unique(catalog_rd_b37$participant_id)
unique_lp_b37 = unique(catalog_rd_b37$platekey)

# catalog b37 is OK, deduplicated, platekey VS participantID
length(unique_pid_b37)
# 99
length(unique_lp_b37)
# 99

# catalog b38
unique_pid_b38 = unique(catalog_rd_b38$participant_id)
unique_lp_b38 = unique(catalog_rd_b38$platekey)

length(unique_pid_b38)
# 76874
length(unique_lp_b38)
# 76881

which_pid_dup_b38 = catalog_rd_b38$participant_id[which(duplicated(catalog_rd_b38$participant_id))]
length(which_pid_dup_b38)
# 7

# select the latest (max) platekey -- associated to the most recent cohort and genome sequenced
dedup_catalog_rd_b38 = catalog_rd_b38 %>%
  filter(!participant_id %in% which_pid_dup_b38)

dup_pid_b38 = catalog_rd_b38 %>% 
  filter(participant_id %in% which_pid_dup_b38) %>%
  group_by(participant_id) %>%
  mutate(latest_lp = max(platekey)) %>%
  ungroup() %>%
  as.data.frame()

dup_pid_b38 = dup_pid_b38 %>% select(latest_lp, participant_id, build, programme)
colnames(dup_pid_b38) = c("platekey", "participant_id", "build", "programme")
dup_pid_b38 = unique(dup_pid_b38)
dim(dup_pid_b38)
# 6  4

dedup_catalog_rd_b38 = rbind(dedup_catalog_rd_b38,
                             dup_pid_b38)
dim(dedup_catalog_rd_b38)
# 76874  4

# check again numbers
length(unique(dedup_catalog_rd_b38$participant_id))
# 76874
length(unique(dedup_catalog_rd_b38$platekey))
# 76874

# Merge catalog_b38 and catalog_b37
rd_catalog = rbind(dedup_catalog_rd_b38,
                   catalog_rd_b37)
dim(rd_catalog)
# 76973  4

length(unique(rd_catalog$platekey))
# 76972
length(unique(rd_catalog$participant_id))
# 76972

# There is one duplicated pid (`50003190`) -- this is in b37 and b38 --> we keep the b38 genome then
to_remove = which(rd_catalog$participant_id %in% "50003190" & rd_catalog$build %in% "GRCh37")
rd_catalog = rd_catalog[-to_remove,]
dim(rd_catalog)
# 76872 4

length(unique(rd_catalog$platekey))
# 76972
length(unique(rd_catalog$participant_id))
# 76972

# Loading last RE clinical data batch (already enriched)
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1124633  28

# Let's only focus on germline in `clin_data`

clin_data = clin_data %>% 
  filter(type %in% c("cancer germline", "experimental germline", "rare disease germline"))
dim(clin_data)
# 1107262  28

# Recover from clin_data genomes not in dedup_cata
l_pid_rd_catalog = unique(rd_catalog$participant_id)

# Just curious, how many do intersect?
length(intersect(clin_data$participant_id, l_pid_rd_catalog))
# 71671

l_in_RE_not_catalog = setdiff(clin_data$participant_id, l_pid_rd_catalog)
length(l_in_RE_not_catalog)
# 15499

# There are 15,499 more participants in RD under V8 RE not in catalog, let's take them
rd_clin_data_not_catalog = clin_data %>%
  filter(participant_id %in% l_in_RE_not_catalog) %>%
  select(platekey, participant_id, genome_build, programme)
dim(rd_clin_data_not_catalog)
# 16251  4

rd_clin_data_not_catalog = unique(rd_clin_data_not_catalog)
dim(rd_clin_data_not_catalog)
# 15779  4

# Before merging, change `rd_catalog$programme` to plural and Uppercase
rd_catalog$programme = rep("Rare Diseases", length(rd_catalog$platekey))
# And change colnames
colnames(rd_clin_data_not_catalog) = colnames(rd_catalog)

rd_catalog_and_RE = rbind(rd_catalog,
                          rd_clin_data_not_catalog)

dim(rd_catalog_and_RE)
# 92751  4

rd_catalog_and_RE = unique(rd_catalog_and_RE)
dim(rd_catalog_and_RE)
# 92751  4

# Are deduplicated?? Double checking again...
length(unique(rd_catalog_and_RE$participant_id))
# 92471
length(unique(rd_catalog_and_RE$platekey))
# 92751

# which duplicated pids??
l_pid_dup = rd_catalog_and_RE[which(duplicated(rd_catalog_and_RE$participant_id)),]$participant_id
length(unique(l_pid_dup))
# 245

# let's filter them out this from rd_catalog_and_RE
dedup_rd_catalog_and_RE = rd_catalog_and_RE %>%
  filter(!participant_id %in% l_pid_dup)
dim(dedup_rd_catalog_and_RE)
# 92226  4

length(unique(dedup_rd_catalog_and_RE$participant_id))
# 92226
length(unique(dedup_rd_catalog_and_RE$platekey))
# 92226

# I've seen that RE V8 is not well enriched with the `type` of the genome
all_germlines = read.csv("~/Documents/STRs/clinical_data/clinical_data/raw/RE_clinical_data_V8/genome_file_paths_and_types_2019-12-04_15-13-29.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(all_germlines)
# 500443  11

# Let's focus only germline genomes
all_germlines = all_germlines %>%
  filter(type %in% c("cancer germline", "experimental germline", "rare disease germline"))

df_pid_dup = clin_data %>% 
  filter(participant_id %in% l_pid_dup) %>%
  select(participant_id, platekey, genome_build, programme)

df_pid_dup = left_join(df_pid_dup,
                       all_germlines %>% 
                         filter(participant_id %in% l_pid_dup) %>% 
                         select(participant_id, platekey, type, genome_build) %>%
                         group_by(participant_id) %>%
                         mutate(latest_platekey = max(platekey)) %>%
                         ungroup() %>%
                         as.data.frame(),
                       by = "participant_id")

# Select only the `cancer germline`, and select the latest ones (there is a mess with b37 and b38...)
df_pid_dup = df_pid_dup %>%
  filter(type %in% "cancer germline") %>%
  select(latest_platekey, participant_id, genome_build.x, programme)

colnames(df_pid_dup) = colnames(dedup_rd_catalog_and_RE)
df_pid_dup = unique(df_pid_dup)
dim(df_pid_dup)
# 242  4

length(unique(df_pid_dup$platekey))
# 242
length(unique(df_pid_dup$participant_id))
# 242


# Include these 242 dedup germline cancer samples to dedup merged table
dedup_rd_catalog_and_RE = rbind(dedup_rd_catalog_and_RE,
                                df_pid_dup)

dim(dedup_rd_catalog_and_RE)
# 92468  4

# Checking again for duplicates...
length(unique(dedup_rd_catalog_and_RE$platekey))
# 92468
length(unique(dedup_rd_catalog_and_RE$participant_id))
# 92468

#------------------------------------------------------------------------------------------------------------------------------
# Let's see if we can "fish" more genomes from earlier batches (EHv2 and EHv3) and population table RE bioinfo group did
# Load popu table
popu_table = read.csv("~/Documents/STRs/clinical_data/clinical_data/aggregate_gvcf_sample_stats_2019-10-03_22-07-31.tsv",
                      sep = "\t",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59356  51

# Is in the population table genomes/participant ids that are missing from our merged Catalog/RE-V8 table?
l_pid_rd_catalog_and_RE = unique(dedup_rd_catalog_and_RE$participant_id)
length(l_pid_rd_catalog_and_RE)
# 92468

length(setdiff(popu_table$participant_id, l_pid_rd_catalog_and_RE))
# 279

# There are 279 PIDs that have been considered when generating population estimations !!  Who are they?
l_pid_popu_not_merged = setdiff(popu_table$participant_id, l_pid_rd_catalog_and_RE)
popu_table_not_merged = popu_table %>%
  filter(participant_id %in% l_pid_popu_not_merged) %>%
  select(platekey, participant_id)
popu_table_not_merged$build = rep("GRCh38", length(popu_table_not_merged$platekey))
popu_table_not_merged$programme = rep("Rare Diseases", length(popu_table_not_merged$platekey))

length(unique(popu_table_not_merged$participant_id))
# 279
length(unique(popu_table_not_merged$platekey))
# 279

dedup_rd_catalog_and_RE = rbind(dedup_rd_catalog_and_RE,
                                popu_table_not_merged)
dim(dedup_rd_catalog_and_RE)
# 92747  4

# Checking again for duplicates...
length(unique(dedup_rd_catalog_and_RE$platekey))
# 92747
length(unique(dedup_rd_catalog_and_RE$participant_id))
# 92747

# EHv2 - summer 2019 batch
l_ehv2_summer = read.table("./batch_august2019_EHv2.5.5/LP_IDs_research_EH_2.5.5.txt", stringsAsFactors = F)
l_ehv2_summer = l_ehv2_summer$V1
length(l_ehv2_summer)
# 86457

# Are all these genomes included in our final dedup merged table?
length(intersect(l_ehv2_summer, dedup_rd_catalog_and_RE$platekey))
# 75730

# which ones are new? and check whether pid is already (in case the are duplicates, with contamination issues etc.)
l_platekeys_in_ehv2_not_merged = setdiff(l_ehv2_summer, dedup_rd_catalog_and_RE$platekey)
length(l_platekeys_in_ehv2_not_merged)
# 10727

# See to which pids they correspond
l_pids_in_ehv2_not_merged = all_germlines %>% 
  filter(platekey %in% l_platekeys_in_ehv2_not_merged) %>%
  select(participant_id) %>%
  unique() %>%
  pull()
length(l_pids_in_ehv2_not_merged)
# 753

new_pids = setdiff(l_pids_in_ehv2_not_merged, dedup_rd_catalog_and_RE$participant_id)
#  111004941 113003384 115011329

# There are only 3 pids that are not included...who are they?
fishing_from_ehv2 = all_germlines %>% 
  filter(participant_id %in% new_pids) %>%
  select(participant_id, platekey, type, genome_build) %>%
  group_by(participant_id) %>%
  mutate(latest_platekey = max(platekey)) %>%
  ungroup() %>%
  as.data.frame()

fishing_from_ehv2 = fishing_from_ehv2 %>%
  select(latest_platekey, participant_id, genome_build)
fishing_from_ehv2 = unique(fishing_from_ehv2)

fishing_from_ehv2$programme = rep("Rare Diseases", length(fishing_from_ehv2$latest_platekey))
colnames(fishing_from_ehv2) = colnames(dedup_rd_catalog_and_RE)


# We recover 3 genomes
dedup_rd_catalog_and_RE = rbind(dedup_rd_catalog_and_RE,
                                fishing_from_ehv2)
dim(dedup_rd_catalog_and_RE)
# 72750  4

# Again checking for duplicates
length(unique(dedup_rd_catalog_and_RE$platekey))
# 92750
length(unique(dedup_rd_catalog_and_RE$participant_id))
# 92750


# EHv3 batch summer / october
l_ehv3_summer = read.table("./batch_october2019_EHv3.1.2/list_platekeys_113911_12oct2019.tsv", stringsAsFactors = F)
l_ehv3_summer = l_ehv3_summer$V1
length(l_ehv3_summer)
#  113911

# Are all these genomes included in our final dedup merged table?
length(intersect(l_ehv3_summer, dedup_rd_catalog_and_RE$platekey))
# 88335

# which ones are new? and check whether pid is already (in case the are duplicates, with contamination issues etc.)
l_platekeys_in_ehv3_not_merged = setdiff(l_ehv3_summer, dedup_rd_catalog_and_RE$platekey)
length(l_platekeys_in_ehv3_not_merged)
# 25576

# See to which pids they correspond
l_pids_in_ehv3_not_merged = all_germlines %>% 
  filter(platekey %in% l_platekeys_in_ehv3_not_merged) %>%
  select(participant_id) %>%
  unique() %>%
  pull()
length(l_pids_in_ehv3_not_merged)
# 860

new_pids = setdiff(l_pids_in_ehv3_not_merged, dedup_rd_catalog_and_RE$participant_id)
#  0
# we already have everything!!

# Some stats before writing final table/genomes into a file
length(unique(dedup_rd_catalog_and_RE$platekey))
# 92750
length(unique(dedup_rd_catalog_and_RE$participant_id))
# 92750

table(dedup_rd_catalog_and_RE$build)
#GRCh37 GRCh38 
#110  77218
table(dedup_rd_catalog_and_RE$programme)
#Cancer Rare Diseases 
#15422         77328 


# Write into a file
# In order to run ExpansionHunter we need the following info for the input file:
# <PLATEKEY>, <PATH TO BAM>, <GENDER>

catalog_rd_b38 = read.csv("./batch_march2020_EHv2.5.5_and_EHv3.2.2/output_catalog_RDb38_280220.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = F)
dim(catalog_rd_b38)
# 76950  8

# RDb37
catalog_rd_b37 = read.csv("./batch_march2020_EHv2.5.5_and_EHv3.2.2/output_catalog_RDb37_280220.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = F)
dim(catalog_rd_b37)
# 99 9

clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1124633  28

l_dedup_rd_catalog_and_RE = unique(dedup_rd_catalog_and_RE$platekey)

# let's retrieve the gender from catalog and clin_data
df_gender = catalog_rd_b38 %>%
  filter(V2 %in% l_dedup_rd_catalog_and_RE) %>%
  select(V2, V5)
colnames(df_gender) = c("platekey", "gender")
df_gender$gender = tolower(df_gender$gender)

df_genderb37 = catalog_rd_b37 %>%
  filter(V2 %in% l_dedup_rd_catalog_and_RE) %>%
  select(V2, V6)
colnames(df_genderb37) = c("platekey", "gender")
df_genderb37$gender = tolower(df_genderb37$gender)

df_gender = rbind(df_gender,
                  df_genderb37)

df_gender_RE = clin_data %>%
  filter(platekey %in% l_dedup_rd_catalog_and_RE) %>%
  select(platekey, participant_phenotypic_sex)
colnames(df_gender_RE) = c("platekey", "gender")

df_gender_RE$gender = gsub("Female", "female", df_gender_RE$gender)
df_gender_RE$gender = gsub("Male", "male", df_gender_RE$gender)


df_gender = rbind(df_gender,
                  df_gender_RE)
df_gender = unique(df_gender)
dim(df_gender)
# 92489

df_gender_popu = popu_table %>%
  filter(platekey %in% l_dedup_rd_catalog_and_RE) %>%
  select(platekey, participant_phenotypic_sex)
colnames(df_gender_popu) = c("platekey", "gender")

df_gender_popu$gender = gsub("Female", "female", df_gender_popu$gender)
df_gender_popu$gender = gsub("Male", "male", df_gender_popu$gender)


df_gender = rbind(df_gender,
                  df_gender_popu)
df_gender = unique(df_gender)
dim(df_gender)
# 92887  2

dedup_rd_catalog_and_RE = left_join(dedup_rd_catalog_and_RE,
                                    df_gender,
                                    by = "platekey")

dedup_rd_catalog_and_RE = unique(dedup_rd_catalog_and_RE)
dim(dedup_rd_catalog_and_RE)
# 92887  5

table(dedup_rd_catalog_and_RE$gender)
#female Indeterminate          male       unknown 
#1         49066            19         43800             1

which(is.na(dedup_rd_catalog_and_RE$gender))
# 0
dedup_rd_catalog_and_RE$platekey[which(dedup_rd_catalog_and_RE$gender == "unknown")]
# "LP3001086-DNA_D10"
dedup_rd_catalog_and_RE$platekey[which(dedup_rd_catalog_and_RE$gender == "")]
# LP3000606-DNA_D12 -- female
dedup_rd_catalog_and_RE$gender[which(dedup_rd_catalog_and_RE$gender == "")] = "female"

l_platekeys_indeterminate = dedup_rd_catalog_and_RE$platekey[which(dedup_rd_catalog_and_RE$gender == "Indeterminate")]
# [1] "LP3001123-DNA_E08" "LP3001008-DNA_H03" "LP3001220-DNA_H07" "LP3001106-DNA_C10" "LP3001163-DNA_F10" "LP3001085-DNA_G04"
#[7] "LP3000840-DNA_D08" "LP3001042-DNA_H12" "LP3001245-DNA_B07" "LP3000828-DNA_E08" "LP3001069-DNA_B08" "LP3001086-DNA_D10"
#[13] "LP3000760-DNA_F09" "LP3001161-DNA_D12" "LP3001248-DNA_E10" "LP3001177-DNA_G01" "LP3001177-DNA_F01" "LP3001152-DNA_F04"
#[19] "LP3001217-DNA_B03"

# these are in catalog b38 and we do have the same info for them, remove all `Indeterminate`
dedup_rd_catalog_and_RE = dedup_rd_catalog_and_RE %>%
  filter(!gender %in% "Indeterminate")
dedup_rd_catalog_and_RE = unique(dedup_rd_catalog_and_RE)
dim(dedup_rd_catalog_and_RE)
# 92867  5

table(dedup_rd_catalog_and_RE$gender)
#female    male unknown 
#49066   43800       1 

# There are 117 genomes/pids with both female and male gender estimations...we need to choose the correct one --> popu_table is taking phenotypic info
# while catalog info is the sex. Let's take the catalog one

l_index_both_gender = which(duplicated(dedup_rd_catalog_and_RE$platekey))
l_platekey_both_gender = dedup_rd_catalog_and_RE$platekey[l_index_both_gender]

df_gender_from_catalog = catalog_rd_b38 %>%
  filter(V2 %in% l_platekey_both_gender) %>%
  select(V2,V5)
dim(df_gender_from_catalog)
# 113  2

# There are 4 missing from there (but existing on the EHv2 summer batch). I do this manually
which(!l_platekey_both_gender %in% df_gender_from_catalog$V2)
# 114 115 116 117
# "LP3000759-DNA_B01" - cancer male"
# LP3000915-DNA_F01" - cancer female
# "LP3000663-DNA_F03" - female 
# "LP3001152-DNA_F09" - male
# Filter out from dedup merged those that are the opposite gender

index_to_remove = which(dedup_rd_catalog_and_RE$platekey %in% "LP3000759-DNA_B01" & dedup_rd_catalog_and_RE$gender == "female")
dedup_rd_catalog_and_RE = dedup_rd_catalog_and_RE[-index_to_remove,]
index_to_remove = which(dedup_rd_catalog_and_RE$platekey %in% "LP3000915-DNA_F01" & dedup_rd_catalog_and_RE$gender == "male")
dedup_rd_catalog_and_RE = dedup_rd_catalog_and_RE[-index_to_remove,]
index_to_remove = which(dedup_rd_catalog_and_RE$platekey %in% "LP3000663-DNA_F03" & dedup_rd_catalog_and_RE$gender == "male")
dedup_rd_catalog_and_RE = dedup_rd_catalog_and_RE[-index_to_remove,]
index_to_remove = which(dedup_rd_catalog_and_RE$platekey %in% "LP3001152-DNA_F09" & dedup_rd_catalog_and_RE$gender == "female")
dedup_rd_catalog_and_RE = dedup_rd_catalog_and_RE[-index_to_remove,]

# change for each of those platekeys
df_gender_from_catalog$V5 = tolower(df_gender_from_catalog$V5)
for (i in 1:length(df_gender_from_catalog$V2)){
  index_to_remove = which(dedup_rd_catalog_and_RE$platekey %in% df_gender_from_catalog$V2[i] & 
                            dedup_rd_catalog_and_RE$gender != df_gender_from_catalog$V5[i])
  dedup_rd_catalog_and_RE = dedup_rd_catalog_and_RE[-index_to_remove,]
  
}

dim(dedup_rd_catalog_and_RE)
# 92750  5

# No duplicates
table(dedup_rd_catalog_and_RE$gender)
# female    male unknown 
# 49000   43749       1 


# I'll define the `unknown` as female..
which(dedup_rd_catalog_and_RE$gender == "unknown")
# 73225
dedup_rd_catalog_and_RE$gender[73225] = "female"
table(dedup_rd_catalog_and_RE$gender)
# female   male 
# 49001  43749

# Let's annotate finally with the BAM path (for the latest platekey)
all_paths_until_feb20 = read.csv("./list_all_genomes_path_together_29022020.tsv",
                                 sep = "\t",
                                 stringsAsFactors = F,
                                 header = F)
dim(all_paths_until_feb20)
# 122586  2

# remove any `plots` or `stats`
all_paths_until_feb20 = all_paths_until_feb20 %>% filter(!grepl("stats", path_bam))
all_paths_until_feb20 = all_paths_until_feb20 %>% filter(!grepl("plots", path_bam))
all_paths_until_feb20 = all_paths_until_feb20 %>% filter(!grepl("time", path_bam))

colnames(all_paths_until_feb20) = c("platekey", "path_bam")

paths_scotland = read.csv("./list_all_genomes_path_together_scotland.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = F)
dim(paths_scotland)
# 5068  2
colnames(paths_scotland) = c("platekey", "path_bam")

all_paths = rbind(all_paths_until_feb20,
                  paths_scotland)

dedup_rd_catalog_and_RE = left_join(dedup_rd_catalog_and_RE,
                                    all_paths,
                                    by = "platekey")
dedup_rd_catalog_and_RE = unique(dedup_rd_catalog_and_RE)

dim(dedup_rd_catalog_and_RE)
# 95929  6

# Take the latest bam_path
dedup_rd_catalog_and_RE = dedup_rd_catalog_and_RE %>%
  group_by(platekey) %>%
  mutate(latest_path = max(path_bam)) %>%
  ungroup() %>%
  as.data.frame()

dedup_rd_catalog_and_RE_full = dedup_rd_catalog_and_RE %>%
  select(platekey, latest_path, gender, build)
dedup_rd_catalog_and_RE_full = unique(dedup_rd_catalog_and_RE_full)
dim(dedup_rd_catalog_and_RE_full)
# 92750  4

# We are missing 1K paths to genomes that do not start with LP (majority from Scotland)
length(which(is.na(dedup_rd_catalog_and_RE_full$latest_path)))
# 81

# write into the file, without the NAs
to_write = dedup_rd_catalog_and_RE_full %>% filter(!is.na(latest_path)) 
dim(to_write)
# 92669  4

# just checking last time
length(to_write$platekey)
# 92669
length(unique(to_write$platekey))
# 92669

# There are many genomes having NAs as build. They must be GRCh38
# AND I assumed WRONGLY that retrieving genomes/data from Catalog RD b38 are GRCh38 build genomes, but no!! they can be also GRCh37!! Augusto asked time ago to merge RD b37 and b38 catalog datasets
# I've got an idea: use the `latest_report` to fix the GRCh38 false ones and change them to GRCh37
upload_report = read.csv("./upload_report.020320.txt",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(upload_report)
# 120648  10

# ACHTUNG!! this table contains everything!! there are duplicated platekeys that have been realigned to GRCh38 a posteriori --> let's take the latest one!!
upload_report = upload_report %>%
  group_by(Platekey) %>%
  mutate(latest_delivery_version = max(Delivery.Version)) %>%
  ungroup() %>%
  select(Platekey, latest_delivery_version) %>%
  as.data.frame()

upload_report = unique(upload_report)
dim(upload_report)
# 115014  2

colnames(upload_report) = c("platekey", "Delivery.Version")

# V1/V2 versions are GRCh37, the rest GRCh38 (as we do have already). Let's focus on GRCh37
to_write = left_join(to_write,
                     upload_report %>% filter(platekey %in% to_write$platekey),
                     by = "platekey")
  
# Recode those `Delivery.version` == V1 or V2 to build=GRCh37
to_write$Delivery.Version = recode(to_write$Delivery.Version, V1 = "GRCh37", V2 = "GRCh37", V4 = "GRCh38")
table(to_write$Delivery.Version)
#  GRCh37  GRCh38 unknown 
# 13024   61835    2276

# what about the `unknown`?
to_write %>% filter(Delivery.Version %in% "unknown") %>% select(build) %>% table()
#GRCh38 
# 1846
to_write$Delivery.Version = recode(to_write$Delivery.Version, unknown = "GRCh38")

# There are still NAs --> GRCh38
to_write$Delivery.Version = to_write$Delivery.Version %>% replace_na("GRCh38")
table(to_write$Delivery.Version)
# GRCh37 GRCh38 
# 13024  79645

to_write_b37 = to_write %>% 
  filter(Delivery.Version %in% "GRCh37") %>% 
  select(platekey, latest_path, gender)
to_write_b37 = unique(to_write_b37)
dim(to_write_b37)
# 13024  3

# GRCh38
to_write_b38 = to_write %>% 
  filter(Delivery.Version %in% "GRCh38") %>% 
  select(platekey, latest_path, gender)
to_write_b38 = unique(to_write_b38)
dim(to_write_b38)
# 79645  3

# Write b37 paths
write.table(to_write_b37, 
            "./batch_march2020_EHv2.5.5_and_EHv3.2.2/list_13024_ouf_of_92669_genomes_GRCh37.csv", 
            sep = ",",
            quote = F, 
            row.names = F,
            col.names = F)

write.table(to_write_b37, 
            "./batch_march2020_EHv2.5.5_and_EHv3.2.2/list_13024_ouf_of_92669_genomes_GRCh37.tsv", 
            sep = "\t",
            quote = F, 
            row.names = F,
            col.names = F)


# Write b38 paths
write.table(to_write_b38, 
            "./batch_march2020_EHv2.5.5_and_EHv3.2.2/list_79645_ouf_of_92669_genomes_GRCh38.csv", 
            sep = ",",
            quote = F, 
            row.names = F,
            col.names = F)

write.table(to_write_b38, 
            "./batch_march2020_EHv2.5.5_and_EHv3.2.2/list_79645_ouf_of_92669_genomes_GRCh38.tsv", 
            sep = "\t",
            quote = F, 
            row.names = F,
            col.names = F)
