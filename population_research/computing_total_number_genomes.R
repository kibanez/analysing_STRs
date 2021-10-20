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

# Let's define as `Neuro` also those having as diseases: Mito or Ultra-rare
# mito and ulutra-rrare should not be included as NOT NEURO cohort
clin_data = clin_data %>% 
  group_by(participant_id) %>%
  mutate(is_neuro = ifelse(rare_diseases_family_id %in% l_fam_neuro | grepl("Mito", diseases_list, ignore.case = T) | grepl("Ultra-rare", diseases_list, ignore.case = T), "Neuro", "NotNeuro")) %>% 
  ungroup() %>%
  as.data.frame()

# Let's check again
# Neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "Neuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 18205 
# Not Neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 36233 

# Breakdown by ancestry for not neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("AFR", "African")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 1239
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("AMR", "American")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 674
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("EUR", "European")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 30754
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("EAS", "East Asian")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 307
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("SAS", "South Asian")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 2765

# Write down into a table to enrich full-mutation and premutation tables (pileup visualisation tables)
write.table(clin_data,
            "~/Documents/STRs/clinical_data/clinical_data/Main_RE_V12_and_Pilot_programmes_enriched_with_neuro_notNeuro.tsv",
            quote = F, col.names = T, row.names = F, sep = "\t")

# Enrich premutation and full-mutation tables with the new neuro - notNeuro info
premut_table = read.csv("./expanded_genomes_main_pilot/feb2021/premutation_table_enriched_with_new_definition_neuro_notNeuro.tsv",
                        header = T, stringsAsFactors = F, sep = "\t")
dim(premut_table)
# 2750  14

patho_table = read.csv("./expanded_genomes_main_pilot/feb2021/full-mutation_table_enriched_with_new_definition_neuro_notNeuro.tsv",
                       stringsAsFactors = F, header = T, sep = "\t")
dim(patho_table)
# 1086  13

premut_table = left_join(premut_table,
                         clin_data %>% select(platekey, is_neuro),
                         by = "platekey")

patho_table = left_join(patho_table,
                        clin_data %>% select(platekey, is_neuro),
                        by = "platekey")

write.table(premut_table, "./expanded_genomes_main_pilot/feb2021/premutation_table_enriched_with_new_definition_neuro_notNeuro2.tsv",
            quote = F, col.names = T, row.names = F, sep = "\t")
write.table(patho_table, "./expanded_genomes_main_pilot/feb2021/full-mutation_table_enriched_with_new_definition_neuro_notNeuro2.tsv",
            quote = F, col.names = T, row.names = F, sep = "\t")
