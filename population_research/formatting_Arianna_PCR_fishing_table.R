# Objective: Arianna worked out all PCR info retrieved from NHNN database
# I need to format the table
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.4"

# Set working directory
setwd("~/Documents/STRs/VALIDATION/PCR_EH_estimations/")

# Load Ari's NHNN PCR table
ari_table = read.csv("~/Documents/STRs/VALIDATION/Arianna_Fishing/STR Repeats in WinPath 06.03.20_per_KRI.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(ari_table)
# 298  10

# load clinical data (in order to match PID with LP_number)
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1124633  28

# merge ari_table with LP_number
ari_table = left_join(ari_table,
                      clin_data %>% select(participant_id, platekey),
                      by = c("PID" = "participant_id"))
ari_table = unique(ari_table)
dim(ari_table)
# 298  11 

# Locus format is a bit diff here
# ATXN1 -- SCA1, SC1B for allele 1 and 2 respectively
# ATXN2 -- SC2A SC2B
# ATXN3 -- SC3A SC3B
# CACNA1A -- SC6A SC6B
# ATXN7 -- SC7A SC7B 
# ATN1 -- ATN1 ATN2
# HTT -- HTT1 HHT2
# C9orf72 -- FTD1 FTD2

table(ari_table$TFC)
#ATN1 ATN2 FTD1 FTD2 HTT1 HTT2 SC1A SC1B SC2A SC2B SC3A SC3B SC6A SC6B SC7A SC7B      
#1    1    7    7    1    1   28   28   28   28   28   28   28   28   28   28 

# First, recode above original locus in the ones we want
# For that we need to have it as factor
ari_table$TFC = as.factor(ari_table$TFC)
ari_table$TFC = recode(ari_table$TFC, 
                       "SC1A" = "ATXN1", "SC1B" = "ATXN1",
                       "SC2A" = "ATXN2", "SC2B" = "ATXN2",
                       "SC3A" = "ATXN3", "SC3B" = "ATXN3",
                       "SC6A" = "CACNA1A", "SC6B" = "CACNA1A",
                       "SC7A" = "ATXN7", "SC7B" = "ATXN7",
                       "ATN1" = "ATN1", "ATN2" = "ATN1",
                       "HTT1" = "HTT", "HTT2" = "HTT",
                       "FTD1" = "C9orf72", "FTD2" = "C9orf72")
                       
# Retrieve list of platekeys together with the locus, in order to take EHv2.5.5 and EHv3.1.2 estimations
df_platekeys_locus= ari_table %>% select(TFC, platekey)
df_platekeys_locus = unique(df_platekeys_locus)
dim(df_platekeys_locus)
# 144 2
write.table(df_platekeys_locus, "../Arianna_Fishing/list_platekeys_locus.tsv", quote = F, row.names = F, col.names = F, sep = "\t")

# Load table Arianna's recoded with PCR_a1 and PCR_a2
final_table_ari = read.csv("../Arianna_Fishing/arianna_table_recoded.tsv",
                           sep = "\t",
                           stringsAsFactors = F,
                           header = T)
dim(final_table_ari)
# 144  6

# Load EHv2 and EHv3 estimations from Arianna's table
eh_final_table = read.csv("../Arianna_Fishing/list_platekeys_locus_EHv2_EHv3.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T)
dim(eh_final_table)
# 144  8

# remove PCR_a1 and PCR_a2
eh_final_table = eh_final_table[,-c(3:4)]

merged_final_table = left_join(final_table_ari,
                               eh_final_table,
                               by = c("platekey" = "platekey", "locus" = "locus"))
dim(merged_final_table)
# 144  10

# Write into a file
write.table(merged_final_table,
            "./NHNN_fishing_Arianna.tsv",
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)
