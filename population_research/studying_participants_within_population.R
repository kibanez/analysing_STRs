# Check which participants we are taking into account in the population research
# are they all unrelated, probands, affected one?
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr)

# set the working directory
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research")

l_participants_pure_ancestry=read.table("list_genomes_56176_pure_ancestry_info.txt",
                                      header = F,
                                      stringsAsFactors = F)
l_participants_pure_ancestry = l_participants_pure_ancestry$V1
length(l_participants_pure_ancestry)
# 56176


re_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_250919.tsv",
                   header = T,
                   stringsAsFactors = F,
                   sep = "\t")
dim(re_data)
# 1056568  26

re_data_subset = re_data %>%
  filter(plate_key.x %in% l_participants_pure_ancestry) %>%
  select(participant_id, plate_key.x, rare_diseases_family_id, participant_type, normalised_specific_disease, genome_build, disease_group, disease_sub_group, specific_disease, programme, family_group_type, affection_status)
dim(re_data_subset)
# 634692  12

length(unique(re_data_subset$plate_key.x))
# 56176

dim(unique(re_data_subset))
# 57704

re_data_subset = unique(re_data_subset)
