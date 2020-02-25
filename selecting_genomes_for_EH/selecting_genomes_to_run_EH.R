# Objective: select genomes (better deduplicated) in order to run afterwards EH through all them
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string # "R version 3.6.2 (2019-12-12)"

# Libraries
library(dplyr)

# Working directory
setwd("~/Documents/STRs/data/research/input")

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
