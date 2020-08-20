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
df_all_genomes = read.csv("~/Documents/STRs/data/research/input/batch_august2020_EHv255_and_EHv322/input/upload_report.2020-08-18.txt",
                          sep = "\t",
                          stringsAsFactors = F,
                          comment.char = "#",
                          header = F)
dim(df_all_genomes)
# 120711  10

# Since HPC data has been moved from Pegasus to Helix, they have changed the name of the platekey
df_all_genomes$V3 = gsub("_copied", "", df_all_genomes$V3)

# Construct new column with the full absolute path to the BAM file
df_all_genomes = df_all_genomes %>%
  mutate(abs_path = paste(paste(paste(V6, "Assembly", sep = "/"), V3, sep = "/"), "bam", sep = "."))

# Select columns
df_all_genomes = df_all_genomes %>%
  select(V3, abs_path, V10)
colnames(df_all_genomes) = c("platekey", "path", "version")

# Create latest or most recent sequenced platekey
df_all_genomes = df_all_genomes %>%
  group_by(platekey) %>%
  mutate(latest_path = max(path)) %>%
  ungroup() %>%
  as.data.frame()

# Load all merged clinical data from RE, MAIN and PILOT together
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_research_cohort/clin_data_merged_V1:V9.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1299442  32

# Let's keep only with the important information and focus on germline genomes
clin_data = clin_data %>%
  filter(type %in% c("cancer germline", "experimental germline", "rare disease", "rare disease germline", "unknown")) %>%
  select(platekey, participant_id, participant_phenotypic_sex, genome_build)
clin_data = unique(clin_data)
dim(clin_data)
# 94365 4

merged_clin_data = clin_data
table(merged_clin_data$genome_build)
#37      38  GRCh37  GRCh38 unknown 

# Put 37 and 38 into GRCh37 and GRCh38
index_37 = which(merged_clin_data$genome_build %in% "37")
index_38 = which(merged_clin_data$genome_build %in% "38")
merged_clin_data$genome_build[index_37] = "GRCh37"
merged_clin_data$genome_build[index_38] = "GRCh38"
table(merged_clin_data$genome_build)
#GRCh37  GRCh38 unknown 

# Load list of genomes/path from March 2020
df_march_b37 = read.csv("./batch_march2020_EHv2.5.5_and_EHv3.2.2/input/list_13024_ouf_of_92669_genomes_GRCh37.csv",
                        stringsAsFactors = F,
                        header = F)
dim(df_march_b37)
# 13024  3

df_march_b38 = read.csv("./batch_march2020_EHv2.5.5_and_EHv3.2.2/input/list_79645_ouf_of_92669_genomes_GRCh38.csv",
                        stringsAsFactors = F,
                        header = F)
dim(df_march_b38)
# 79645  3

l_platekeys_march = c(df_march_b37$V1,
                      df_march_b38$V1)
length(l_platekeys_march)
# 92669

l_pid_march = merged_clin_data %>%
  filter(platekey %in% l_platekeys_march) %>%
  select(participant_id) %>%
  unique() %>%
  pull()
length(l_pid_march)
# 88359

# Enrich with build info
df_march_b37$build = rep("GRCh37", length(df_march_b37$V1))
df_march_b38$build = rep("GRCh38", length(df_march_b38$V1))

# The ones we know are the latest ones, and then from PIDs, we can get the PLATEKEY and the GENDER
# <PLATEKEY>,<PATH>,<GENDER>,<BUILD>
df_final_list = rbind(df_march_b37,
                      df_march_b38)
dim(df_final_list)
# 92669  4

colnames(df_final_list) = c("platekey", "path", "gender", "build")

l_final_platekey = unique(df_final_list$platekey)
length(l_final_platekey)
# 92669

# Load total unique genomes included in batch1 and batch2 of population analysis
l_popu_genomes = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/list_79849_unique_genomes_batch1_batch2.txt",
                            stringsAsFactors = F)
l_popu_genomes = l_popu_genomes$V1
length(l_popu_genomes)
# 79849

# Check which ones are NEW to what we had in batch march 2020
l_popu_genomes_new = setdiff(l_popu_genomes,
                             l_final_platekey)
length(l_popu_genomes_new)
# 404

l_pid_popu = merged_clin_data %>%
  filter(platekey %in% l_popu_genomes_new) %>%
  select(participant_id) %>%
  unique() %>%
  pull()
length(l_pid_popu)
# 399

# Deduplicate  l_pid_popu
df_popu = merged_clin_data %>%
  filter(platekey %in% l_popu_genomes_new)
df_popu = unique(df_popu)
dim(df_popu)
# 470  4

which_pid_popu_dup = df_popu$participant_id[which(duplicated(df_popu$participant_id))]
length(which_pid_popu_dup)
# 71

# Many of them is because platekey has been sequenced in b37 and b38, let's see if removing genome_build we still have dups
df_popu = df_popu %>%
  select(platekey, participant_id, participant_phenotypic_sex)
df_popu = unique(df_popu)
dim(df_popu)
# 400  3

which_pid_popu_dup = df_popu$participant_id[which(duplicated(df_popu$participant_id))]
length(which_pid_popu_dup)
# 1

df_popu %>% filter(participant_id %in% which_pid_popu_dup)
#platekey participant_id participant_phenotypic_sex
#1 LP3000651-DNA_E02      215001136                       <NA>
#2 LP3000651-DNA_E02      215001136                     Female

# Let's remove the NA one
index_to_remove = which(df_popu$participant_id %in% which_pid_popu_dup & is.na(df_popu$participant_phenotypic_sex))
df_popu = df_popu[-index_to_remove,]

# Check again for duplicates
which_pid_popu_dup = df_popu$participant_id[which(duplicated(df_popu$participant_id))]
length(which_pid_popu_dup)
# 0

# Number of PIDs and Platekeys
length(unique(df_popu$platekey))
# 399
length(unique(df_popu$participant_id))
# 399

# Enrich df_popu with path
df_popu = left_join(df_popu,
                    df_all_genomes %>% select(platekey, latest_path),
                    by = "platekey")
df_popu = unique(df_popu)
dim(df_popu)
# 399 4

# Reorder columns
df_popu = df_popu %>% select(platekey, latest_path, participant_phenotypic_sex)
df_popu$build = rep("GRCh38", length(df_popu$platekey))
colnames(df_popu) = colnames(df_final_list)

# Update final list of genomes - updating with the NEW ONES, just double checking this again
new_platekeys = setdiff(df_popu$platekey,
                        df_final_list$platekey)
length(new_platekeys)
# 399

df_final_list = rbind(df_final_list,
                      df_popu)

df_final_list$gender = tolower(df_final_list$gender)
dim(df_final_list)
# 93068  4
df_final_list = unique(df_final_list)
dim(df_final_list)
# 93068  4

length(unique(df_final_list$platekey))
# 93068

# Check whether the list of genomes in popu are included
list_extra_popu_genomes = read.table("~/Documents/STRs/data/research/input/list_644_genomes_not_included_in_batch_march2020_but_popu.tsv",
                                     stringsAsFactors = F,
                                     header = F)
list_extra_popu_genomes = list_extra_popu_genomes$V1
length(list_extra_popu_genomes)
# 644

length(intersect(list_extra_popu_genomes, 
                 df_final_list$platekey))
# 259

which_new = setdiff(list_extra_popu_genomes,
                    df_final_list$platekey)
length(which_new)
# 385

# Are these somatic??
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_research_cohort/clin_data_merged_V1:V9.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1299442  32

clin_data %>% filter(platekey %in% which_new) %>% select(type) %>% table()
# Pilot

# They are 385 new platekeys that have no clinical information, let's assign them "female" even though they are not
df_new_no_gender = df_all_genomes %>% filter(platekey %in% which_new)
length(unique(df_new_no_gender$platekey))
# 385

# Recode version to GRCh37/GRCh38
df_new_no_gender = df_new_no_gender %>%
  mutate(build = case_when(version == "V2" ~ "GRCh37",
                           version == "V1" ~ "GRCh37",
                           version == "V4" ~ "GRCh38"))

# Assign female to all, even though is not true - but we don't have info
df_new_no_gender$gender = rep("female", length(df_new_no_gender$platekey))
# Order 
df_new_no_gender = df_new_no_gender %>%
  select(platekey, latest_path, gender, build)
df_new_no_gender = unique(df_new_no_gender)
dim(df_new_no_gender)
# 391  4

# Remove duplicates of having b37 and b38 genome builds
l_dups = df_new_no_gender$platekey[which(duplicated(df_new_no_gender$platekey))]
length(l_dups)
# 6

final_df_new_no_gender = df_new_no_gender %>%
  filter(!platekey %in% l_dups)

df_new_no_gender = df_new_no_gender %>%
  filter(platekey %in% l_dups)
length(unique(df_new_no_gender$platekey))
# 6

# They are in both b37 and b38 --> let's take b38
df_new_no_gender = df_new_no_gender %>%
  filter(build %in% "GRCh38")
length(unique(df_new_no_gender$platekey))
# 6

# merge them after dedup
final_df_new_no_gender = rbind(final_df_new_no_gender,
                               df_new_no_gender)
final_df_new_no_gender = unique(final_df_new_no_gender)
dim(final_df_new_no_gender)
# 385  4

length(unique(final_df_new_no_gender$platekey))
# 385

# Reorder columns
colnames(final_df_new_no_gender) = colnames(df_final_list)

# Merge all
df_final_list = rbind(df_final_list,
                      final_df_new_no_gender)
df_final_list = unique(df_final_list)
dim(df_final_list)
# 93453  4

length(unique(df_final_list$platekey))
# 93453

to_write_b37 = df_final_list %>% 
  filter(build %in% "GRCh37") %>% 
  select(platekey, path, gender)
to_write_b37 = unique(to_write_b37)
dim(to_write_b37)
# 13401  3

# GRCh38
to_write_b38 = df_final_list %>% 
  filter(build %in% "GRCh38") %>% 
  select(platekey, path, gender)
to_write_b38 = unique(to_write_b38)
dim(to_write_b38)
# 80052  3

# Write b37 paths
write.table(to_write_b37, 
            "./batch_august2020_EHv255_and_EHv322/list_13401_ouf_of_93453_genomes_GRCh37.csv", 
            sep = ",",
            quote = F, 
            row.names = F,
            col.names = F)

write.table(to_write_b37, 
            "./batch_august2020_EHv255_and_EHv322/list_13401_ouf_of_93453_genomes_GRCh37.tsv", 
            sep = "\t",
            quote = F, 
            row.names = F,
            col.names = F)


# Write b38 paths
write.table(to_write_b38, 
            "./batch_august2020_EHv255_and_EHv322/list_80052_ouf_of_93453_genomes_GRCh38.csv", 
            sep = ",",
            quote = F, 
            row.names = F,
            col.names = F)

write.table(to_write_b38, 
            "./batch_august2020_EHv255_and_EHv322/list_80052_ouf_of_93453_genomes_GRCh38.tsv", 
            sep = "\t",
            quote = F, 
            row.names = F,
            col.names = F)
