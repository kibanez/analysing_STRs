# Objective: define genomes to include in the validation of long-read sequencing: pacbio and nanopore
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv", sep = "\t", header = T, stringsAsFactors = F)
dim(clin_data)
# 1124633  28

# Load the list of genomes proposed from Neuroscience group (6th floor)
list_long = read.table("~/Documents/STRs/ANALYSIS/long_read/Validation_cohort/list_participantID.txt", stringsAsFactors = F)
list_long = list_long$V1
length(list_long)
# 28

# there is 1 duplicate
length(unique(list_long))
# 27

# Subset of 28 genomes in the clinical data and take the important information
clin_data_subset = clin_data %>% 
  filter(participant_id %in% list_long) %>% 
  select(participant_id, platekey, file_path, type, rare_diseases_family_id, plate_key, biological_relationship_to_proband, 
         participant_type, normalised_specific_disease, genome_build, genetic_vs_reported_results, disease_group, disease_sub_group, 
         specific_disease, participant_phenotypic_sex, programme_consent_status, programme, family_group_type)

dim(clin_data_subset)
# 1163  18

clin_data_subset = unique(clin_data_subset)
dim(clin_data_subset)
# 27  18