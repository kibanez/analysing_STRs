# Objective: select genomes (better deduplicated) in order to run afterwards EH through all them
# The idea is to merge or join all data available to us: RE batch + population aggregated gVCF (in case we have more genomes to fish) + earlier EH batches
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
all_genomes = read.table("~/Documents/STRs/data/research/input/list_all_genomes_25022020.tsv", stringsAsFactors = F)
all_genomes = all_genomes$V1
length(all_genomes)
# 121506

# we do this once, and after this, we will load the whole table with genomes and path to the genomes
df_all_genomes = c()
for(i in 1:length(all_genomes)){
  path_to_genome = all_genomes[i]
  bam_name = strsplit(path_to_genome, "/")[[1]][length(strsplit(path_to_genome, "/")[[1]])-2]
  df_all_genomes = rbind(df_all_genomes,
                         cbind(bam_name, path_to_genome))
  
}

df_all_genomes = as.data.frame(df_all_genomes)
write.table(df_all_genomes,
            "~/Documents/STRs/data/research/input/list_all_genomes_path_together_25022020.tsv",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

# Loading last RE clinical data batch (already enriched)
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1124633  28

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
