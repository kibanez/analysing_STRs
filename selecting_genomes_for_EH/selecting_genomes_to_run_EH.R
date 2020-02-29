# Objective: select genomes (better deduplicated) in order to run afterwards EH through all them
# The idea is to merge or join all data available to us: RE batch + population aggregated gVCF (in case we have more genomes to fish) + earlier EH batches
# RE batch (V8 last batch from November)
# Population aggregated gVCF file
# EHv2 batch from summer 2019
# EHv3 batch from summer 2019

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

# Loading all genomes sequenced by GEL so far (25th Feb 2020)
#all_genomes = read.table("~/Documents/STRs/data/research/input/list_all_genomes_25022020.tsv", stringsAsFactors = F)
#all_genomes = all_genomes$V1
#length(all_genomes)
# 121506

# we do this once, and after this, we will load the whole table with genomes and path to the genomes
#df_all_genomes = c()
#for(i in 1:length(all_genomes)){
#  path_to_genome = all_genomes[i]
#  bam_name = strsplit(path_to_genome, "/")[[1]][length(strsplit(path_to_genome, "/")[[1]])-2]
#  df_all_genomes = rbind(df_all_genomes,
#                         cbind(bam_name, path_to_genome))
#}

#df_all_genomes = as.data.frame(df_all_genomes)
#write.table(df_all_genomes,
#            "~/Documents/STRs/data/research/input/list_all_genomes_path_together_25022020.tsv",
#            quote = F,
#            col.names = F,
#            row.names = F,
#            sep = "\t")

# Load table with all genomes together with their paths
df_all_genomes = read.csv("~/Documents/STRs/data/research/input/list_all_genomes_path_together_25022020.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = F)
dim(df_all_genomes)
# 121506  2

# Platekey, PID info retrieved from Catalog
# RD b38
catalog_rd_b38 = read.csv("./batch_march2020_EHv2.5.5_and_EHv3.1.2/output_catalog_RDb38_280220.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = F)
dim(catalog_rd_b38)
# 76950  8

# RDb37
catalog_rd_b37 = read.csv("./batch_march2020_EHv2.5.5_and_EHv3.1.2/output_catalog_RDb37_280220.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = F)
dim(catalog_rd_b37)
# 99 9

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

# Recover from clin_data genomes not in dedup_cata
l_rd_catalog = 

# There are duplicates here, more genomes than participantIDs
length(unique(clin_data$platekey))
# 104843

length(unique(clin_data$participant_id))
# 87395

all_germlines = read.csv("~/Documents/STRs/clinical_data/clinical_data/raw/RE_clinical_data_V8/genome_file_paths_and_types_2019-12-04_15-13-29.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(all_germlines)
# 500443  11

table(all_germlines$type)
#cancer germline        cancer somatic experimental germline  experimental somatic rare disease germline 
#73813                 52249                   268                   159                373954 


l_germline_pid = all_germlines %>%
  filter(type %in% c("cancer germline", "rare disease germline")) %>%
  select(participant_id) %>%
  unique() %>%
  pull()
length(l_germline_pid)
# 87168

# Load popu table
popu_table = read.csv("~/Documents/STRs/clinical_data/clinical_data/aggregate_gvcf_sample_stats_2019-10-03_22-07-31.tsv",
                      sep = "\t",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59356  51

popu_platekeys = unique(popu_table$platekey)
length(popu_platekeys)
# 59356
length(intersect(popu_platekeys, all_germlines$platekey))
# 58441

setdiff(popu_platekeys, all_germlines$platekey)
