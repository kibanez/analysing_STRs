# Objective: study and retrieve numbers re demographic within 100kGP
# unrelated genomes, filtering out 125bp genomes
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(tidyverse)

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/")

# Load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V10_and_Pilot_programmes.tsv",
                     stringsAsFactors = F, 
                     header = T,
                     sep = "\t")
dim(clin_data)
# 2101385  24

# Load unrel genomes
l_unrel = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/l_unrelated_55603_genomes_batch2.txt",
                     stringsAsFactors = F)
l_unrel = l_unrel$V1
length(l_unrel)
# 55603

# Load list of 2x125bp sequenced genomes
l_125 = read.table("./list_genomes_125bp_100kGP.tsv", stringsAsFactors = F)
l_125 = l_125$V1
length(l_125)
# 15830

# List of unrel genomes 2x150bp
l_unrel_no125 = setdiff(l_unrel, l_125)
length(l_unrel_no125)
# 54437

clin_data_paper = clin_data %>%
  filter(platekey %in% l_unrel_no125) %>%
  select(platekey, participant_phenotypic_sex, year_of_birth, participant_ethnic_category,
         programme, family_group_type, affection_status, superpopu)
clin_data_paper = unique(clin_data_paper)
dim(clin_data_paper)
# 54437  8

clin_data_paper = clin_data_paper %>%
  group_by(platekey) %>%
  mutate(age = 2020 - as.integer(year_of_birth)) %>%
  ungroup() %>%
  as.data.frame()
