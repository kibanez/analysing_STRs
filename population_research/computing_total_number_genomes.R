# Script to compute total number of genomes in the 100kGP
# We filter out genomes sequenced in read-length 125bp
# All vs not considering neuro
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/")

l_125 = read.table("./list_genomes_125bp_100kGP.tsv", stringsAsFactors = F)
l_125 = l_125$V1

l_unrel = read.table("./l_unrelated_55603_genomes_batch2.txt", stringsAsFactors = F)
l_unrel = l_unrel$V1

# Total number of genomes, filtering out 125bp sequenced genomes
l_unrel_not125 = unique(setdiff(l_unrel, l_125))
length(l_unrel_not125)
# 54437

# Total number of genomes excluding Neuro
# Update October 2021: we want to consider as Neuro, patients that have been assigned "Mito" or "ultra-rare" diseases

# Load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V12_and_Pilot_programmes.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(clin_data)
# 2472865  26

# First, compute of 54437 how many they are RD and Cancer
clin_data %>% filter(platekey %in% l_unrel_not125) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 54437
clin_data %>% filter(platekey %in% l_unrel_not125, programme %in% "Cancer") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 14628 (26.87%)
clin_data %>% filter(platekey %in% l_unrel_not125, programme %in% "Rare Diseases") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 39809 (73.13%)

l_fam_neuro = clin_data %>%
  filter(grepl("neuro", diseasegroup_list, ignore.case = T)) %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull()
length(l_fam_neuro)
# 14717

clin_data = clin_data %>% 
  filter(platekey %in% l_unrel_not125) %>%
  group_by(rare_diseases_family_id) %>% 
  mutate(is_neuro = ifelse(rare_diseases_family_id %in% l_fam_neuro, "Neuro", "NotNeuro")) %>% 
  ungroup() %>% 
  as.data.frame() 

# Neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "Neuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 17608
# Not Neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 36830

# Let's define as `NotNeuro` also those having as diseases: Mito or Ultra-rare
clin_data = clin_data %>% 
  group_by(participant_id) %>%
  filter(grepl("Mito", diseases_list, ignore.case = T) | grepl("Ultra-rare", diseases_list, ignore.case = T)) %>%
  mutate(is_neuro = "Neuro") %>%
  ungroup() %>%
  as.data.frame()

# Let's check again
# Neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "Neuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 610 
# Not Neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 0 


clin_data_notNeuro = clin_data %>%
  filter(!grepl("neuro", diseasegroup_list, ignore.case = TRUE))

clin_data_notNeuro %>% filter(platekey %in% l_unrel_not125) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 48656


clin_data_notNeuro = clin_data %>%
  filter(!grepl("neuro", diseasegroup_list, ignore.case = TRUE) | !grepl("Mito", diseases_list,ignore.case = TRUE) | !grepl("Ultra-rare", diseases_list,ignore.case = TRUE))

  select(platekey) %>%
  unique() %>%
  pull()