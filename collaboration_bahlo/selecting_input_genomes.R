# Script to generate the list of Unrelated PROBANDS not neurological
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/FAME3/")

l_unrel = read.table("./l_unrelated_55603_genomes_batch2.txt", stringsAsFactors = F)
l_unrel = l_unrel$V1

# Load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V12_and_Pilot_programmes.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(clin_data)
# 2472865  26

# Filter out Neurological people (probands and their familials)
# Including `Mito` and `ultra-rare` under neuro
l_fam_neuro = clin_data %>%
  filter(grepl("neuro", diseasegroup_list, ignore.case = T)) %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull()
length(l_fam_neuro)
# 14717

l_fam_neuro = c(l_fam_neuro,
                clin_data %>%
                  filter(grepl("Mito", diseases_list, ignore.case = T)) %>%
                  select(rare_diseases_family_id) %>%
                  unique() %>%
                  pull())

l_fam_neuro = c(l_fam_neuro,
                clin_data %>%
                  filter(grepl("Ultra-rare", diseases_list, ignore.case = T)) %>%
                  select(rare_diseases_family_id) %>%
                  unique() %>%
                  pull())
length(l_fam_neuro)
# 16492

clin_data = clin_data %>% 
  filter(platekey %in% l_unrel) %>%
  group_by(rare_diseases_family_id) %>% 
  mutate(is_neuro = ifelse(rare_diseases_family_id %in% l_fam_neuro, "Neuro", "NotNeuro")) %>% 
  ungroup() %>% 
  as.data.frame() 

# Neuro
clin_data %>% filter(platekey %in% l_unrel, is_neuro %in% "Neuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 20243
# Not Neuro
clin_data %>% filter(platekey %in% l_unrel, is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 35361

# We want to have Not Neuro genomes that are PROBANDS
l_not_neuro = clin_data %>% filter(platekey %in% l_unrel, is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull()
length(l_not_neuro)
# 35361

l_not_neuro_proband = clin_data %>%
  filter(platekey %in% l_not_neuro, participant_type %in% "Proband") %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_not_neuro_proband)
# 9404
