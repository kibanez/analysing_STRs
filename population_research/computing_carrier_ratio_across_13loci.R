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

# Load the whole table for 100kGP - case-controls 
table_100cc_QC = read.csv("./table_platekey_locus_QC_inspection.tsv",
                          stringsAsFactors = F,
                          header = T,
                          sep = "\t")
dim(table_100cc_QC)
# 1878  13

# count unique PID included in the cases_controls (i.e. what is the total number of genomes that we have data on)
total_number_of_participants_analysed <- length(unique(clin_data$participant_id))
# 88826

# For each locus, add a new column to `clin_data` if the repeat size of each locus is larger than path threshold
l_locus = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "C9ORF72", "DMPK", "FMR1", "FXN", "HTT", "TBP")
l_patho_cutoff = c(38,48,44,33,60,36,60,20,50,66,49)

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


