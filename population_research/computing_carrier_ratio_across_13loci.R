# Objective: from the work we have done by inspecting visually all pileups, compute the carrier ratio for each locus
# unrelated, probands unrelated, probands not neuro
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/")

# Load clinical data for 90,863 genomes for which we do have EHv3 estimations across 13 loci
clin_data = read.csv("~/Documents/STRs/ANALYSIS/cases_controls/batch_march/EHv322/table_13_loci_across_90863_genomes_each_row_allele_EHv3.2.2_90863platekeys_88827pids.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(clin_data)
# 91857  46

# remove platekeys that have no PIDs
clin_data  = clin_data %>% filter(participant_id != ".")
dim(clin_data)
# 89821  46

# Load list of platekeys from GenQA
l_platekeys_genQA = read.table("/Users/kibanez/Documents/STRs/VALIDATION/genQA/genQA/list_platekeys_b1_total_RD_probands_and_cancer_expanded_after_QC_locus_merged_genQA.txt",
                             stringsAsFactors = F,
                             header = F)
l_platekeys_genQA = unique(l_platekeys_genQA$V1)
length(l_platekeys_genQA)
# 51

# Include new column called `genQA` as boolean
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(genQA = ifelse(platekey %in% l_platekeys_genQA, TRUE, FALSE))
table(clin_data$genQA)
# ALL FALSE, Good, we don't want to include them as part of 100K cohort

# Load the whole table for 100kGP - case-controls 
table_100cc_QC = read.csv("./table_platekey_locus_QC_inspection.tsv",
                          stringsAsFactors = F,
                          header = T,
                          sep = "\t")
dim(table_100cc_QC)
# 1783  13

# count unique PID included in the cases_controls (i.e. what is the total number of genomes that we have data on)
total_number_of_participants_analysed <- length(unique(clin_data$participant_id))
# 88826

# For each locus, add a new column to `clin_data` if the repeat size of each locus is larger than path threshold
#l_locus = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "C9ORF72", "DMPK", "FMR1", "FXN", "HTT", "TBP")
l_locus = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "C9ORF72", "DMPK", "FXN","TBP")
#l_patho_cutoff = c(38,48,44,33,60,36,60,20,50,66,49)

# AR
# First thing, for FMR1 and AR, we need to transform NA values in a1|a2 to 0
values_AR_a1 = which(is.na(clin_data$AR_a1))
# 0
values_AR_a2 = which(is.na(clin_data$AR_a2))
clin_data$AR_a2[values_AR_a2] = 0

clin_data = clin_data %>% 
  mutate(AR_before_VI = if_else(AR_a1 >= 38 | AR_a2 >= 38, TRUE, FALSE))

# ATN1
clin_data = clin_data %>% 
  mutate(ATN1_before_VI = if_else(ATN1_a1 >= 48 | ATN1_a2 >= 48, TRUE, FALSE))

# ATXN1
clin_data = clin_data %>% 
  mutate(ATXN1_before_VI = if_else(ATXN1_a1 >= 44 | ATXN1_a2 >= 44, TRUE, FALSE))

# ATXN2
clin_data = clin_data %>% 
  mutate(ATXN2_before_VI = if_else(ATXN2_a1 >= 33 | ATXN2_a2 >= 33, TRUE, FALSE))

# ATXN3
clin_data = clin_data %>% 
  mutate(ATXN3_before_VI = if_else(ATXN3_a1 >= 60 | ATXN3_a2 >= 60, TRUE, FALSE))

# ATXN7
clin_data = clin_data %>% 
  mutate(ATXN7_before_VI = if_else(ATXN7_a1 >= 36 | ATXN7_a2 >= 36, TRUE, FALSE))

# CACNA1A
clin_data = clin_data %>% 
  mutate(CACNA1A_before_VI = if_else(CACNA1A_a1 >= 20 | CACNA1A_a2 >= 20, TRUE, FALSE))

# C9ORF72
clin_data = clin_data %>% 
  mutate(C9ORF72_before_VI = if_else(C9ORF72_a1 >= 60 | C9ORF72_a2 >= 60, TRUE, FALSE))

# DMPK
clin_data = clin_data %>% 
  mutate(DMPK_before_VI = if_else(DMPK_a1 >= 50 | DMPK_a2 >= 50, TRUE, FALSE))

# HTT
# Arianna's table

# FMR1
# There are no genomes with expansions beyond pathogenic cutoff

# FXN
clin_data = clin_data %>% 
  mutate(FXN_before_VI = if_else(FXN_a1 >= 66 | FXN_a2 >= 66, TRUE, FALSE))

# TBP
clin_data = clin_data %>% 
  mutate(TBP_before_VI = if_else(TBP_a1 >= 49 | TBP_a2 >= 49, TRUE, FALSE))


# Now let's annotate or enrich with after visual inspection data
# AR
l_platekeys_AR_true_after = table_100cc_QC %>% filter(locus %in% "AR", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_AR_true_after)
# 69

# ATN1
l_platekeys_ATN1_true_after = table_100cc_QC %>% filter(locus %in% "ATN1", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_ATN1_true_after)
# 8

# ATXN1
l_platekeys_ATXN1_true_after = table_100cc_QC %>% filter(locus %in% "ATXN1", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_ATXN1_true_after)
# 17

# ATXN2
l_platekeys_ATXN2_true_after = table_100cc_QC %>% filter(locus %in% "ATXN2", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_ATXN2_true_after)
# 38

# ATXN3
l_platekeys_ATXN3_true_after = table_100cc_QC %>% filter(locus %in% "ATXN3", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_ATXN3_true_after)
# 1

# ATXN7
l_platekeys_ATXN7_true_after = table_100cc_QC %>% filter(locus %in% "ATXN7", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_ATXN7_true_after)
# 10 

# CACNA1A
l_platekeys_CACNA1A_true_after = table_100cc_QC %>% filter(locus %in% "CACNA1A", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_CACNA1A_true_after)
# 23

# C9ORF72
l_platekeys_C9ORF72_true_after = table_100cc_QC %>% filter(locus %in% "C9ORF72", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_C9ORF72_true_after)
# 26 

# DMPK
l_platekeys_DMPK_true_after = table_100cc_QC %>% filter(locus %in% "DMPK", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_DMPK_true_after)
# 63

# FXN
l_platekeys_FXN_true_after = table_100cc_QC %>% filter(locus %in% "FXN", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_FXN_true_after)
# 1165 (1 genQA or more - check)

# TBP
l_platekeys_TBP_true_after = table_100cc_QC %>% filter(locus %in% "TBP", Final.decision %in% "Yes") %>% select(platekey) %>% unique() %>% pull()
length(l_platekeys_TBP_true_after)
# 2 (1 in GenQA)

# Include new columns after visual inspection
# AR
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(AR_after_VI = ifelse(platekey %in% l_platekeys_AR_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(AR_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 64

# ATN1
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(ATN1_after_VI = ifelse(platekey %in% l_platekeys_ATN1_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(ATN1_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 5

# ATXN1
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(ATXN1_after_VI = ifelse(platekey %in% l_platekeys_ATXN1_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(ATXN1_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 17

# ATXN2
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(ATXN2_after_VI = ifelse(platekey %in% l_platekeys_ATXN2_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(ATXN2_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 38

# ATXN3
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(ATXN3_after_VI = ifelse(platekey %in% l_platekeys_ATXN3_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(ATXN3_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 1

# ATXN7
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(ATXN7_after_VI = ifelse(platekey %in% l_platekeys_ATXN7_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(ATXN7_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 9

# CACNA1A
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(CACNA1A_after_VI = ifelse(platekey %in% l_platekeys_CACNA1A_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(CACNA1A_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 22

# C9ORF72
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(C9ORF72_after_VI = ifelse(platekey %in% l_platekeys_C9ORF72_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(C9ORF72_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 25

# DMPK
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(DMPK_after_VI = ifelse(platekey %in% l_platekeys_DMPK_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(DMPK_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 59

# FXN
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(FXN_after_VI = ifelse(platekey %in% l_platekeys_FXN_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(FXN_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 1120

# TBP
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(TBP_after_VI = ifelse(platekey %in% l_platekeys_TBP_true_after, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()
clin_data %>% filter(TBP_after_VI) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 1

dim(clin_data)
# 89821  68

# 1 - PROBANDS in RD or Cancer
#select genomes from the cases_control table that are only probands in Rare disease or cancer (=N/A in main, Proband in Pilot, all cancer germline)
clin_data_RD_probands_and_cancer = clin_data %>% 
  filter(biological_relationship_to_proband %in% "N/A" | 
           biological_relationship_to_proband %in% "Proband" | 
           is.na(biological_relationship_to_proband) | 
           programme %in% "Cancer")
dim(clin_data_RD_probands_and_cancer)
# 49774  69

#count unique PID included in the probands rd and cancer table (i.e. how many people are RD or cancer probands)
total_number_of_pids_RD_probands_and_cancer_analysed = length(unique(clin_data_RD_probands_and_cancer$participant_id))
# 49447

# Compute carrier ratio for each locus
# select ONLY genomes that have an expansion that passed visual QC in the RD and cancer probands
# Create a df for PROBANDS (cancer and rd)
df_probands = data.frame()
for (i in 1:length(l_locus)){
  locus_after_VI = paste(l_locus[i], "after_VI", sep = "_")
  total_RD_probands_and_cancer_expanded_after_QC_locus = clin_data_RD_probands_and_cancer %>% 
    filter(eval(parse(text=locus_after_VI)) == TRUE) %>%
    select(participant_id) %>%
    unique() %>%
    pull() %>%
    length()
  
  freq_carrier = round(total_number_of_pids_RD_probands_and_cancer_analysed / total_RD_probands_and_cancer_expanded_after_QC_locus,digits = 2)
  ratio_freq_carrier = paste("1 in", as.character(freq_carrier), sep = " ")
  
  ci_max = round(total_number_of_pids_RD_probands_and_cancer_analysed/(total_number_of_pids_RD_probands_and_cancer_analysed*((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed)-1.96*sqrt((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed)*(1-total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed)/total_number_of_pids_RD_probands_and_cancer_analysed))), digits = 2)
  ci_min = round(total_number_of_pids_RD_probands_and_cancer_analysed/(total_number_of_pids_RD_probands_and_cancer_analysed*((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed)+1.96*sqrt((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed)*(1-total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed)/total_number_of_pids_RD_probands_and_cancer_analysed))), digits = 2)
  
  ci_ratio= as.character(paste(as.character(ci_min), as.character(ci_max), sep = "-"))
  
  df_probands = rbind(df_probands,
                      cbind(l_locus[i], total_RD_probands_and_cancer_expanded_after_QC_locus, total_number_of_pids_RD_probands_and_cancer_analysed, ratio_freq_carrier, ci_ratio))
}
# write into a table
write.table(df_probands,
            "total_numbers_and_ratio_PROBANDS_in_RD_or_Cancer.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# Now the same, but focusing on EUR genomes
total_number_of_pids_RD_probands_and_cancer_analysed_EUR = clin_data_RD_probands_and_cancer %>% filter(superpopu %in% "EUR") %>% select(participant_id) %>% unique() %>% pull() %>% length()
df_probands_eur = data.frame()
for (i in 1:length(l_locus)){
  locus_after_VI = paste(l_locus[i], "after_VI", sep = "_")
  total_RD_probands_and_cancer_expanded_after_QC_locus = clin_data_RD_probands_and_cancer %>% 
    filter(eval(parse(text=locus_after_VI)) == TRUE, superpopu %in% "EUR") %>%
    select(participant_id) %>%
    unique() %>%
    pull() %>%
    length()
  
  freq_carrier = round(total_number_of_pids_RD_probands_and_cancer_analysed_EUR / total_RD_probands_and_cancer_expanded_after_QC_locus,digits = 2)
  ratio_freq_carrier = paste("1 in", as.character(freq_carrier), sep = " ")
  
  ci_max = round(total_number_of_pids_RD_probands_and_cancer_analysed_EUR/(total_number_of_pids_RD_probands_and_cancer_analysed_EUR*((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed_EUR)-1.96*sqrt((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed_EUR)*(1-total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed_EUR)/total_number_of_pids_RD_probands_and_cancer_analysed_EUR))), digits = 2)
  ci_min = round(total_number_of_pids_RD_probands_and_cancer_analysed_EUR/(total_number_of_pids_RD_probands_and_cancer_analysed_EUR*((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed_EUR)+1.96*sqrt((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed_EUR)*(1-total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_analysed_EUR)/total_number_of_pids_RD_probands_and_cancer_analysed_EUR))), digits = 2)
  
  ci_ratio= as.character(paste(as.character(ci_min), as.character(ci_max), sep = "-"))
  
  df_probands_eur = rbind(df_probands_eur,
                      cbind(l_locus[i], total_RD_probands_and_cancer_expanded_after_QC_locus, total_number_of_pids_RD_probands_and_cancer_analysed_EUR, ratio_freq_carrier, ci_ratio))
}
# write into a table
write.table(df_probands_eur,
            "total_numbers_and_ratio_PROBANDS_in_RD_or_Cancer_EUR.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# 2 - PROBANDS in RD or Cancer - BUT NOT NEURO!!!
# Select those genomes/participants that have NOT been recruited in Neurological disorders
clin_data_RD_probands_and_cancer_notNeuro = clin_data_RD_probands_and_cancer %>% 
  filter(!grepl("neuro", diseasegroup_list, ignore.case = TRUE))
dim(clin_data_RD_probands_and_cancer_notNeuro)
# 35636  69

total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed = length(unique(clin_data_RD_probands_and_cancer_notNeuro$participant_id))
# 35450
 
# Compute carrier ratio for each locus
# select ONLY genomes that have an expansion that passed visual QC in the RD and cancer probands
# Create a df for PROBANDS_NOT_NEURO (cancer and rd)
df_probands_notNeuro = data.frame()
for (i in 1:length(l_locus)){
  locus_after_VI = paste(l_locus[i], "after_VI", sep = "_")
  total_RD_probands_and_cancer_expanded_after_QC_locus = clin_data_RD_probands_and_cancer_notNeuro %>% 
    filter(eval(parse(text=locus_after_VI)) == TRUE) %>%
    select(participant_id) %>%
    unique() %>%
    pull() %>%
    length()
  
  freq_carrier = round(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed / total_RD_probands_and_cancer_expanded_after_QC_locus,digits = 2)
  ratio_freq_carrier = paste("1 in", as.character(freq_carrier), sep = " ")
  
  ci_max = round(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed/(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed*((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed)-1.96*sqrt((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed)*(1-total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed)/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed))), digits = 2)
  ci_min = round(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed/(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed*((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed)+1.96*sqrt((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed)*(1-total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed)/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed))), digits = 2)
  
  ci_ratio= as.character(paste(as.character(ci_min), as.character(ci_max), sep = "-"))
  
  df_probands_notNeuro = rbind(df_probands_notNeuro,
                               cbind(l_locus[i], total_RD_probands_and_cancer_expanded_after_QC_locus, total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed, ratio_freq_carrier, ci_ratio))
}
# write into a table
write.table(df_probands_notNeuro,
            "total_numbers_and_ratio_PROBANDS_in_RD_notNeuro_or_Cancer.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# Now the same, but focusing on EUR genomes
total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR = clin_data_RD_probands_and_cancer_notNeuro %>% filter(superpopu %in% "EUR") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 27068
df_probands_notNeuro_eur = data.frame()
for (i in 1:length(l_locus)){
  locus_after_VI = paste(l_locus[i], "after_VI", sep = "_")
  total_RD_probands_and_cancer_expanded_after_QC_locus = clin_data_RD_probands_and_cancer_notNeuro %>% 
    filter(eval(parse(text=locus_after_VI)) == TRUE,
           superpopu %in% "EUR") %>%
    select(participant_id) %>%
    unique() %>%
    pull() %>%
    length()
  
  freq_carrier = round(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR / total_RD_probands_and_cancer_expanded_after_QC_locus,digits = 2)
  ratio_freq_carrier = paste("1 in", as.character(freq_carrier), sep = " ")
  
  ci_max = round(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR/(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR*((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR)-1.96*sqrt((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR)*(1-total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR)/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR))), digits = 2)
  ci_min = round(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR/(total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR*((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR)+1.96*sqrt((total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR)*(1-total_RD_probands_and_cancer_expanded_after_QC_locus/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR)/total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR))), digits = 2)
  
  ci_ratio= as.character(paste(as.character(ci_min), as.character(ci_max), sep = "-"))
  
  df_probands_notNeuro_eur = rbind(df_probands_notNeuro_eur,
                                   cbind(l_locus[i], total_RD_probands_and_cancer_expanded_after_QC_locus, total_number_of_pids_RD_probands_and_cancer_notNeuro_analysed_EUR, ratio_freq_carrier, ci_ratio))
}

# write into a table
write.table(df_probands_notNeuro_eur,
            "total_numbers_and_ratio_PROBANDS_in_RD_notNeuro_or_Cancer_EUR.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

