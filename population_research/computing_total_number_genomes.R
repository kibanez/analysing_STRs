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

# Dec'21: we spotted out that we also need to consider as Neuro families with the proband recruited under ultra-rare or mito
# Let's define as `Neuro` also those having as diseases: Mito or Ultra-rare
# mito and ulutra-rrare should not be included as NOT NEURO cohort
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
  filter(platekey %in% l_unrel_not125) %>%
  group_by(rare_diseases_family_id) %>% 
  mutate(is_neuro = ifelse(rare_diseases_family_id %in% l_fam_neuro, "Neuro", "NotNeuro")) %>% 
  ungroup() %>% 
  as.data.frame() 

# Neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "Neuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 19771 (before 17608)
# Not Neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 34667 (before 36830)

# Check also by gender (for AR)
# Not Neuro
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Female") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 18970 (before: 19800)
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Male") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 15696 (before: 16431)

# Breakdown by popu
# ALL
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Female") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 29151
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Male") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 25286

# EUR
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Female", superpopu %in% "EUR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 24695
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Male", superpopu %in% "EUR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 21085

# AFR
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Female", superpopu %in% "AFR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 926
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Male", superpopu %in% "AFR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 839

# AMR
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Female", superpopu %in% "AMR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 572
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Male", superpopu %in% "AMR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 578

# EAS
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Female", superpopu %in% "EAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 257
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Male", superpopu %in% "EAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 160

# SAS
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Female", superpopu %in% "SAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 2323
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, participant_phenotypic_sex %in% "Male", superpopu %in% "SAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 2266

# NOT NEURO
# ALL
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Female") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 18970 (before: 19800)
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Male") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 15696 (before: 16432)

# EUR
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Female", superpopu %in% "EUR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 16259 (before: 16926)
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Male", superpopu %in% "EUR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 13237 (before: 13827)

# AFR
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Female", superpopu %in% "AFR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 642 (before: 655)
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Male", superpopu %in% "AFR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 569 (before: 584)

# AMR
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Female", superpopu %in% "AMR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 316 (before:335)
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Male", superpopu %in% "AMR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 320 (before: 339)

# EAS
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Female", superpopu %in% "EAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 186 (before: 194)
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Male", superpopu %in% "EAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 108 (before: 113)

# SAS
# XX
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Female", superpopu %in% "SAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 1318 (before: 1432)
# XY
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", participant_phenotypic_sex %in% "Male", superpopu %in% "SAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 1234 (before: 1333)


# Breakdown by ancestry for not neuro
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("AFR", "African")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 1211 (before: 1239)
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("AMR", "American")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 636 (before: 674)
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("EUR", "European")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 29497 (before: 30754)
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("EAS", "East Asian")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 294 (before: 307)
clin_data %>% filter(platekey %in% l_unrel_not125, is_neuro %in% "NotNeuro", superpopu %in% c("SAS", "South Asian")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 2552 (before: 2765)

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
                         clin_data %>% select(platekey, is_neuro, participant_phenotypic_sex, year_of_birth, programme),
                         by = "platekey")
premut_table = unique(premut_table)
patho_table = left_join(patho_table,
                        clin_data %>% select(platekey, is_neuro, participant_phenotypic_sex, year_of_birth, programme),
                        by = "platekey")
patho_table = unique(patho_table)

write.table(premut_table, "./expanded_genomes_main_pilot/feb2021/premutation_table_enriched_with_new_definition_neuro_notNeuro2_with_gender_yob.tsv",
            quote = F, col.names = T, row.names = F, sep = "\t")
write.table(patho_table, "./expanded_genomes_main_pilot/feb2021/full-mutation_table_enriched_with_new_definition_neuro_notNeuro2_with_gender_yob.tsv",
            quote = F, col.names = T, row.names = F, sep = "\t")

# For DMPK we changed or updated the premutation cut-off to >43
premut_table_DMPK = read.csv("./expanded_genomes_main_pilot/feb2021/premutation_table_DMPK_enriched_with_new_definition_neuro_notNeuro.tsv",
                             stringsAsFactors = F, header = T, sep = "\t")
dim(premut_table_DMPK)
# 272 10
premut_table_DMPK = left_join(premut_table_DMPK,
                              clin_data %>% select(platekey, is_neuro, participant_phenotypic_sex, year_of_birth, programme),
                              by = "platekey")
premut_table_DMPK = unique(premut_table_DMPK)
dim(premut_table_DMPK)
# 272  12

write.table(premut_table_DMPK, "./expanded_genomes_main_pilot/feb2021/premutation_table_DMPK_enriched_with_new_definition_neuro_notNeuro2_with_gender.tsv",
            quote = F, col.names = T, row.names = F, sep = "\t")
