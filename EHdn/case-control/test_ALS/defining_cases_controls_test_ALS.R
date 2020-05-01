# Objective: define and track here the modus operandi when selecting samples for cases and controls for PILOT adult ataxia
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.9.0/case-control/analysis/test_ALS/")

# Load list of cases
l_cases = read.table("./list_ALS_GRCh38.txt", stringsAsFactors = F)
l_cases = l_cases$V1


# Load Main data
main_clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_300320.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T)
dim(main_clin_data)
# 1124633  31

# Cases - `l_cases`
l_main_cases = l_cases

# Controls
# - ONLY probands
# - rare disease germline genomes (exclude for now cancer genomes)
# - `disease_group` NOT in neurology
# - for `main` take GRCh38
main_controls = main_clin_data %>%
  filter(participant_type %in% "Proband",
         programme %in% "Rare Diseases",
         genome_build %in% "GRCh38",
         year_of_birth < 2000,
         !plate_key %in% l_main_cases,
         !grepl("[Nn][Ee][Uu][Rr][Oo]", main_clin_data$disease_group))
dim(main_controls)
# 124282  31

length(unique(main_controls$plate_key))
# 11355
